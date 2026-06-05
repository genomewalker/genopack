#include "run_ingest.hpp"
#include <genopack/archive.hpp>
#include <genopack/xqual_columns.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/section_checksum.hpp>
#include <genopack/toc.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace genopack {

namespace {

std::vector<std::filesystem::path> collect_gpk_paths(const std::filesystem::path& pack_path) {
    std::vector<std::filesystem::path> paths;
    if (pack_path.extension() == ".gpk") { paths.push_back(pack_path); return paths; }
    for (const auto& e : std::filesystem::directory_iterator(pack_path))
        if (e.path().extension() == ".gpk") paths.push_back(e.path());
    std::sort(paths.begin(), paths.end());
    return paths;
}

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> out;
    size_t start = 0;
    while (true) {
        size_t tab = line.find('\t', start);
        if (tab == std::string::npos) { out.push_back(line.substr(start)); break; }
        out.push_back(line.substr(start, tab - start));
        start = tab + 1;
    }
    for (auto& s : out) if (!s.empty() && s.back() == '\r') s.pop_back();
    return out;
}

// Strip a FASTA-ish extension so a CheckM2 "Name" (basename of the input FASTA)
// resolves to the archive accession: GCF_x.fna(.gz) -> GCF_x.
std::string strip_fasta_ext(std::string s) {
    static const char* exts[] = {".gz", ".fasta", ".fna", ".faa", ".fa"};
    bool changed = true;
    while (changed) {
        changed = false;
        for (const char* e : exts) {
            const size_t el = std::string(e).size();
            if (s.size() > el && s.compare(s.size() - el, el, e) == 0) {
                s.erase(s.size() - el); changed = true;
            }
        }
    }
    return s;
}

// Case-insensitive header match against any of the candidate names.
int find_col(const std::vector<std::string>& hdr, std::initializer_list<const char*> names) {
    auto lc = [](std::string x) { for (auto& c : x) c = static_cast<char>(::tolower(c)); return x; };
    for (size_t i = 0; i < hdr.size(); ++i) {
        std::string h = lc(hdr[i]);
        for (const char* n : names) if (h == lc(n)) return static_cast<int>(i);
    }
    return -1;
}

bool parse_float(const std::string& s, float& out) {
    if (s.empty() || s == "NA" || s == "NA\r" || s == "None" || s == "-") return false;
    try { out = std::stof(s); return true; } catch (...) { return false; }
}

// Accumulate a (axis,tool,method) measurement for one genome into the columns map,
// keyed by column identity so a re-ingest overwrites the same cell.
struct ColAcc {
    ColumnKey                           key;
    std::unordered_map<uint64_t, float> values;
};

ColAcc& column_for(std::unordered_map<uint64_t, ColAcc>& cols, const ColumnKey& key) {
    const uint64_t id = column_identity_hash(key);
    auto it = cols.find(id);
    if (it == cols.end()) it = cols.emplace(id, ColAcc{key, {}}).first;
    return it->second;
}

// Resolve a TSV row name to a genome_id in this archive (exact, then ext-stripped).
const GenomeId* resolve(const std::unordered_map<std::string, GenomeId>& exact,
                        const std::unordered_map<std::string, GenomeId>& stripped,
                        const std::string& name) {
    auto it = exact.find(name);
    if (it != exact.end()) return &it->second;
    it = stripped.find(strip_fasta_ext(name));
    if (it != stripped.end()) return &it->second;
    return nullptr;
}

// Append a fresh SEC_XQAL (dropping any prior one) to a single .gpk.
void write_xqual(const std::filesystem::path& gpk, const std::vector<XqualColumn>& columns) {
    int lock_fd = ::open(gpk.c_str(), O_RDWR);
    if (lock_fd < 0) throw std::runtime_error("ingest: cannot open for locking: " + gpk.string());
    struct LockGuard { int fd; explicit LockGuard(int f):fd(f){::flock(fd,LOCK_EX);}
                       ~LockGuard(){::flock(fd,LOCK_UN);::close(fd);} } lock(lock_fd);

    MmapFileReader mmap;
    mmap.open(gpk);
    auto toc = TocReader::read(mmap);
    const auto* tail = mmap.ptr_at<TailLocator>(mmap.size() - sizeof(TailLocator));
    uint64_t prev_toc_offset = tail->toc_offset;
    uint64_t generation      = toc.header.generation + 1;
    mmap.close();

    AppendWriter writer;
    writer.open_append(gpk);
    uint64_t section_id = toc.next_section_id();
    SectionDesc xqual_sd = xqual_write(writer, section_id, columns);

    TocWriter new_toc;
    for (const auto& sd : toc.sections)
        if (sd.type != SEC_XQAL) new_toc.add_section(sd);
    new_toc.add_section(xqual_sd);

    writer.flush();
    stamp_section_checksums(gpk.c_str(), new_toc.sections(), /*only_if_zero=*/true);

    new_toc.finalize(writer, generation,
                     toc.header.live_genome_count, toc.header.total_genome_count,
                     prev_toc_offset, toc.header.catalog_root_section_id,
                     toc.header.accession_root_section_id, toc.header.tombstone_root_section_id);
    spdlog::info("ingest: wrote XQAL section_id={} ({} columns) to {}",
                 xqual_sd.section_id, columns.size(), gpk.filename().string());
}

// Parse one external tool TSV, adding measurements into `cols` for genomes present
// in this archive. Returns rows matched.
size_t ingest_tsv(const std::filesystem::path& tsv,
                  const char* tool,
                  const char* method,
                  std::initializer_list<const char*> name_cols,
                  std::initializer_list<const char*> compl_cols,
                  std::initializer_list<const char*> contam_cols,
                  const std::unordered_map<std::string, GenomeId>& exact,
                  const std::unordered_map<std::string, GenomeId>& stripped,
                  std::unordered_map<uint64_t, ColAcc>& cols)
{
    std::ifstream f(tsv);
    if (!f) throw std::runtime_error("ingest: cannot open " + tsv.string());
    std::string line;
    if (!std::getline(f, line)) return 0;
    auto hdr = split_tab(line);
    const int ni = find_col(hdr, name_cols);
    const int ci = find_col(hdr, compl_cols);
    const int ki = find_col(hdr, contam_cols);
    if (ni < 0) throw std::runtime_error("ingest: no genome-name column in " + tsv.string());
    if (ci < 0 && ki < 0)
        throw std::runtime_error("ingest: no completeness/contamination column in " + tsv.string());

    ColumnKey compl_key{qual_axis::COMPLETENESS, tool, method, Unit::Fraction01, 1, 0, 0, 0,
                        std::string(tool) + " completeness"};
    ColumnKey contam_key{qual_axis::CONTAMINATION, tool, method, Unit::Fraction01, 1, 0, 0, 0,
                        std::string(tool) + " contamination"};

    size_t matched = 0;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        auto f_ = split_tab(line);
        if (static_cast<int>(f_.size()) <= ni) continue;
        const GenomeId* gid = resolve(exact, stripped, f_[ni]);
        if (!gid) continue;
        ++matched;
        float v;
        if (ci >= 0 && static_cast<int>(f_.size()) > ci && parse_float(f_[ci], v))
            column_for(cols, compl_key).values[*gid] = v / 100.0f;
        if (ki >= 0 && static_cast<int>(f_.size()) > ki && parse_float(f_[ki], v))
            column_for(cols, contam_key).values[*gid] = v / 100.0f;
    }
    return matched;
}

} // namespace

int cmd_ingest(const std::filesystem::path& pack_path,
               const std::filesystem::path& checkm2_tsv,
               const std::filesystem::path& anvio_tsv)
{
    if (checkm2_tsv.empty() && anvio_tsv.empty())
        throw std::runtime_error("ingest: provide at least one of --checkm2 / --anvio");

    auto gpks = collect_gpk_paths(pack_path);
    if (gpks.empty())
        throw std::runtime_error("ingest: no .gpk found at " + pack_path.string());

    for (const auto& gp : gpks) {
        // Build name -> genome_id maps (exact + extension-stripped).
        std::unordered_map<std::string, GenomeId> exact, stripped;
        std::unordered_map<uint64_t, ColAcc> cols;
        {
            ArchiveReader ar;
            ar.open(gp);
            ar.scan_genome_accessions([&](std::string_view acc, GenomeId gid) {
                if (ar.is_deleted(gid)) return;
                exact.emplace(std::string(acc), gid);
                stripped.emplace(strip_fasta_ext(std::string(acc)), gid);
            });
            // Seed from any existing XQAL so a second ingest accumulates.
            if (ar.has_xqual()) {
                for (auto& xc : xqual_read_all(*ar.xqual_reader())) {
                    ColAcc acc{std::move(xc.key), std::move(xc.values)};
                    cols.emplace(column_identity_hash(acc.key), std::move(acc));
                }
            }
            ar.close();
        }
        if (exact.empty()) { spdlog::info("ingest: {} empty, skipping", gp.filename().string()); continue; }

        size_t m_cm2 = 0, m_anv = 0;
        if (!checkm2_tsv.empty())
            m_cm2 = ingest_tsv(checkm2_tsv, qual_tool::CHECKM2, qual_method::ML,
                               {"Name", "genome", "accession", "bin"},
                               {"Completeness", "completeness"},
                               {"Contamination", "contamination"},
                               exact, stripped, cols);
        if (!anvio_tsv.empty())
            m_anv = ingest_tsv(anvio_tsv, qual_tool::ANVIO, qual_method::MARKER_SCG,
                               {"bin name", "bin_name", "genome", "name", "accession"},
                               {"percent_completion", "% completion", "completion", "completeness"},
                               {"percent_redundancy", "% redundancy", "redundancy", "contamination"},
                               exact, stripped, cols);

        if (cols.empty()) {
            spdlog::info("ingest: {} no matching genomes, skipping", gp.filename().string());
            continue;
        }

        std::vector<XqualColumn> columns;
        columns.reserve(cols.size());
        for (auto& [_, acc] : cols)
            columns.push_back(XqualColumn{std::move(acc.key), std::move(acc.values)});

        write_xqual(gp, columns);
        spdlog::info("ingest: {} — CheckM2 matched {}, anvi'o matched {}, {} columns total",
                     gp.filename().string(), m_cm2, m_anv, columns.size());
    }
    return 0;
}

} // namespace genopack
