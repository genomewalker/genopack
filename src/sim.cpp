#include <genopack/sim.hpp>
#include <genopack/util.hpp>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <numeric>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace genopack {

namespace {

// Fixed-size chunking (CheckM2-compatible).
std::vector<std::string>
chunk_contigs(const std::vector<std::pair<std::string, std::string>>& contigs,
              int chunk_size, int min_chunk)
{
    std::vector<std::string> chunks;
    const auto cs = static_cast<size_t>(chunk_size);
    const auto mc = static_cast<size_t>(min_chunk);
    for (const auto& [name, seq] : contigs)
        for (size_t pos = 0; pos < seq.size(); pos += cs) {
            size_t len = std::min(cs, seq.size() - pos);
            if (len >= mc) chunks.push_back(seq.substr(pos, len));
        }
    return chunks;
}

// Lognormal chunking: each chunk length is sampled from lognormal(ln(n50), sigma).
// Mimics real MAG assemblies where N50 ~ n50 and contig lengths span 2–3 orders of magnitude.
// Redundancy from any taxon sharing markers will appear proportional to shared SCG content.
std::vector<std::string>
chunk_contigs_lognormal(const std::vector<std::pair<std::string, std::string>>& contigs,
                        int n50, double sigma, int min_chunk, std::mt19937_64& rng)
{
    std::vector<std::string> chunks;
    const double mu = std::log(static_cast<double>(n50));
    const size_t mc = static_cast<size_t>(min_chunk);
    std::normal_distribution<double> nd(mu, sigma);
    for (const auto& [name, seq] : contigs) {
        size_t pos = 0;
        while (pos < seq.size()) {
            size_t cs = static_cast<size_t>(std::max(
                static_cast<double>(min_chunk),
                std::exp(nd(rng))));
            cs = std::min(cs, seq.size() - pos);
            if (cs >= mc) chunks.push_back(seq.substr(pos, cs));
            pos += cs;
        }
    }
    return chunks;
}

uint64_t total_bp(const std::vector<std::string>& chunks)
{
    uint64_t bp = 0;
    for (const auto& c : chunks) bp += c.size();
    return bp;
}

std::string base_accession(const std::filesystem::path& p)
{
    std::string s = p.filename().string();
    for (const char* ext : {".fa.gz", ".fna.gz", ".fasta.gz", ".gz",
                             ".fasta", ".fna", ".fa"}) {
        size_t el = std::string(ext).size();
        if (s.size() >= el && s.compare(s.size() - el, el, ext) == 0)
            s.resize(s.size() - el);
    }
    return s;
}

void write_fasta(const std::filesystem::path& path,
                 const std::vector<std::pair<std::string, std::string>>& contigs)
{
    std::ofstream out(path);
    if (!out) throw std::runtime_error("sim: cannot write " + path.string());
    std::string buf;
    buf.reserve(1 << 20);
    for (const auto& [name, seq] : contigs) {
        buf += '>'; buf += name; buf += '\n';
        for (size_t i = 0; i < seq.size(); i += 80) {
            buf.append(seq, i, std::min<size_t>(80, seq.size() - i));
            buf += '\n';
        }
    }
    out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
}

// A job is one output file: (ref, comp, cont_level, contam_source, rep).
// contam_idx == -1 means no contamination (cont_frac must be 0.0).
struct Job { int ref_idx, comp_idx, cont_idx, contam_idx, rep; };

struct GtRow {
    std::string synth_acc, source_acc, taxonomy;
    std::string contam_source_acc, contam_label, contam_taxonomy;
    std::filesystem::path fa_path;
    int      rep;
    double   comp_target, comp_actual, cont_target, cont_actual;
    uint32_t n_contigs;
    uint64_t total_bp_out;
};

} // namespace

std::vector<std::pair<std::string, std::string>>
parse_fasta_named(const std::string& fasta_buf)
{
    std::vector<std::pair<std::string, std::string>> out;
    std::string name, seq;
    auto flush = [&] {
        if (!name.empty()) out.emplace_back(std::move(name), std::move(seq));
        name.clear(); seq.clear();
    };
    for (size_t i = 0; i < fasta_buf.size(); ) {
        if (fasta_buf[i] == '>') {
            flush();
            size_t j = i + 1;
            while (j < fasta_buf.size() && fasta_buf[j] != '\n' &&
                   fasta_buf[j] != ' ' && fasta_buf[j] != '\t' && fasta_buf[j] != '\r')
                ++j;
            name.assign(fasta_buf, i + 1, j - (i + 1));
            size_t eol = fasta_buf.find('\n', i);
            i = (eol == std::string::npos) ? fasta_buf.size() : eol + 1;
        } else {
            size_t eol = fasta_buf.find('\n', i);
            if (eol == std::string::npos) eol = fasta_buf.size();
            for (size_t k = i; k < eol; ++k)
                if (fasta_buf[k] != '\r') seq += fasta_buf[k];
            i = eol + 1;
        }
    }
    flush();
    return out;
}

std::vector<std::pair<std::string, std::string>>
make_fragment(const std::vector<std::pair<std::string, std::string>>& contigs,
              double comp_frac, int chunk_size, int min_chunk,
              int contig_n50, double contig_sigma,
              std::mt19937_64& rng)
{
    std::vector<std::string> chunks = (contig_n50 > 0)
        ? chunk_contigs_lognormal(contigs, contig_n50, contig_sigma, min_chunk, rng)
        : chunk_contigs(contigs, chunk_size, min_chunk);
    std::vector<std::string> kept;
    if (comp_frac >= 1.0) {
        kept = std::move(chunks);
    } else {
        const uint64_t tot    = total_bp(chunks);
        const uint64_t target = static_cast<uint64_t>(
            std::llround(static_cast<double>(tot) * comp_frac));
        std::shuffle(chunks.begin(), chunks.end(), rng);
        uint64_t acc = 0;
        for (auto& c : chunks) {
            if (acc >= target) break;
            acc += c.size();
            kept.push_back(std::move(c));
        }
    }
    std::vector<std::pair<std::string, std::string>> frag;
    frag.reserve(kept.size());
    for (size_t i = 0; i < kept.size(); ++i)
        frag.emplace_back("ctg_" + std::to_string(i), std::move(kept[i]));
    return frag;
}

double
mix_contamination(std::vector<std::pair<std::string, std::string>>& fragments,
                  const std::vector<std::pair<std::string, std::string>>& contam_chunks,
                  double cont_frac, std::mt19937_64& rng)
{
    if (cont_frac <= 0.0 || contam_chunks.empty()) return 0.0;
    uint64_t host_bp = 0;
    for (const auto& [n, s] : fragments) host_bp += s.size();
    const double denom        = std::max(1e-9, 1.0 - cont_frac);
    const uint64_t tgt_contam = static_cast<uint64_t>(
        std::llround(static_cast<double>(host_bp) * cont_frac / denom));
    std::vector<size_t> order(contam_chunks.size());
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), rng);
    const size_t base = fragments.size();
    uint64_t added = 0;
    for (size_t k = 0; k < order.size() && added < tgt_contam; ++k) {
        const auto& seq = contam_chunks[order[k]].second;
        fragments.emplace_back("ctg_" + std::to_string(base + k), seq);
        added += seq.size();
    }
    return static_cast<double>(added) / static_cast<double>(host_bp + added);
}

int run_sim(const SimConfig& cfg)
{
    namespace fs = std::filesystem;
    fs::create_directories(cfg.output_dir);
    const fs::path gt_path  = cfg.output_tsv.empty()
        ? cfg.output_dir / "sim_manifest.tsv" : cfg.output_tsv;
    const fs::path add_path = cfg.manifest_tsv.empty()
        ? cfg.output_dir / "add_manifest.tsv" : cfg.manifest_tsv;

    // Load and chunk all reference genomes
    std::vector<std::vector<std::pair<std::string, std::string>>> ref_contigs(cfg.refs.size());
    std::vector<std::string> ref_acc(cfg.refs.size());
    std::vector<uint64_t>    ref_total_bp(cfg.refs.size(), 0);
    for (size_t i = 0; i < cfg.refs.size(); ++i) {
        ref_contigs[i] = parse_fasta_named(decompress_gz(cfg.refs[i].fasta));
        ref_acc[i]     = base_accession(cfg.refs[i].fasta);
        if (ref_contigs[i].empty())
            throw std::runtime_error("sim: empty reference " + cfg.refs[i].fasta.string());
        for (const auto& [n, s] : ref_contigs[i]) ref_total_bp[i] += s.size();
    }

    // Load and chunk all contamination sources
    // contam_pools[k] = pre-chunked (renamed cx_N) entries for contam source k
    const int n_contams = static_cast<int>(cfg.contams.size());
    std::vector<std::vector<std::pair<std::string, std::string>>> contam_pools(n_contams);
    std::vector<std::string> contam_acc(n_contams);
    for (int k = 0; k < n_contams; ++k) {
        auto raw = chunk_contigs(
            parse_fasta_named(decompress_gz(cfg.contams[k].fasta)),
            cfg.chunk_size, cfg.min_chunk);
        contam_pools[k].reserve(raw.size());
        for (size_t i = 0; i < raw.size(); ++i)
            contam_pools[k].emplace_back("cx_" + std::to_string(i), std::move(raw[i]));
        contam_acc[k] = base_accession(cfg.contams[k].fasta);
    }

    // Build job list.
    // At cont_frac=0: one job per (ref, comp, rep) regardless of contam sources.
    // At cont_frac>0 with no contam sources: skip (can't contaminate without a source).
    // At cont_frac>0 with contam sources: one job per (ref, comp, cont_level, contam_src, rep).
    std::vector<Job> jobs;
    for (int ri = 0; ri < (int)cfg.refs.size(); ++ri)
        for (int ci = 0; ci < (int)cfg.completeness.size(); ++ci)
            for (int xi = 0; xi < (int)cfg.contamination.size(); ++xi) {
                const double cont = cfg.contamination[xi];
                if (cont <= 0.0) {
                    for (int rep = 0; rep < cfg.reps; ++rep)
                        jobs.push_back({ri, ci, xi, -1, rep});
                } else {
                    for (int ki = 0; ki < std::max(1, n_contams); ++ki)
                        for (int rep = 0; rep < cfg.reps; ++rep)
                            if (ki < n_contams)  // skip if no contam sources
                                jobs.push_back({ri, ci, xi, ki, rep});
                }
            }

    std::vector<GtRow> rows(jobs.size());
    spdlog::info("sim: {} refs × {} comp × {} cont × {} contam_srcs × {} reps = {} jobs",
                 cfg.refs.size(), cfg.completeness.size(),
                 cfg.contamination.size(), std::max(1, n_contams), cfg.reps, jobs.size());

#ifdef _OPENMP
    omp_set_num_threads(std::max(1, cfg.threads));
    #pragma omp parallel for schedule(dynamic)
#endif
    for (long j = 0; j < (long)jobs.size(); ++j) {
        const Job& job    = jobs[j];
        const double comp = cfg.completeness[job.comp_idx];
        const double cont = cfg.contamination[job.cont_idx];
        const int    ki   = job.contam_idx;

        // Deterministic per-job seed
        uint64_t s = cfg.seed;
        s = s * 1000003u + (uint64_t)(job.ref_idx    + 1);
        s = s * 1000003u + (uint64_t)(job.comp_idx   + 1);
        s = s * 1000003u + (uint64_t)(job.cont_idx   + 1);
        s = s * 1000003u + (uint64_t)((ki < 0 ? 0 : ki) + 1);
        s = s * 1000003u + (uint64_t)(job.rep         + 1);
        std::mt19937_64 rng(s);

        auto frag = make_fragment(ref_contigs[job.ref_idx], comp,
                                  cfg.chunk_size, cfg.min_chunk,
                                  cfg.contig_n50, cfg.contig_sigma, rng);
        uint64_t host_bp = 0;
        for (const auto& [n, sq] : frag) host_bp += sq.size();
        const double comp_actual = static_cast<double>(host_bp)
                                 / static_cast<double>(ref_total_bp[job.ref_idx]);

        const double cont_actual = (ki >= 0)
            ? mix_contamination(frag, contam_pools[ki], cont, rng)
            : 0.0;

        uint64_t out_bp = 0;
        for (const auto& [n, sq] : frag) out_bp += sq.size();

        // Filename: sim_REF_r00_c050_x05[_LABEL|_k00].fa
        // No contam suffix when cont=0 (clean genome, source irrelevant)
        std::string contam_tag;
        if (ki >= 0) {
            const std::string& lbl = cfg.contams[ki].label;
            contam_tag = "_" + (lbl.empty() ? "k" + std::to_string(ki) : lbl);
        }

        char fname[320];
        std::snprintf(fname, sizeof(fname), "sim_%s_r%02d_c%03d_x%02d%s.fa",
                      ref_acc[job.ref_idx].c_str(), job.rep,
                      (int)std::lround(comp * 100.0),
                      (int)std::lround(cont * 100.0),
                      contam_tag.c_str());

        const fs::path fa = cfg.output_dir / fname;
        write_fasta(fa, frag);

        std::string synth(fname, std::strlen(fname) - 3);
        GtRow& r           = rows[j];
        r.synth_acc        = synth;
        r.source_acc       = ref_acc[job.ref_idx];
        r.taxonomy         = cfg.refs[job.ref_idx].taxonomy;
        r.contam_source_acc= ki >= 0 ? contam_acc[ki] : "";
        r.contam_label     = ki >= 0 ? cfg.contams[ki].label : "";
        r.contam_taxonomy  = ki >= 0 ? cfg.contams[ki].taxonomy : "";
        r.fa_path          = fa;
        r.rep              = job.rep;
        r.comp_target      = comp;
        r.comp_actual      = comp_actual;
        r.cont_target      = cont;
        r.cont_actual      = cont_actual;
        r.n_contigs        = static_cast<uint32_t>(frag.size());
        r.total_bp_out     = out_bp;
    }

    {
        std::ofstream tsv(gt_path);
        if (!tsv) throw std::runtime_error("sim: cannot write " + gt_path.string());
        tsv << "synth_acc\tsource_acc\trep\tcomp_target\tcomp_actual\t"
               "cont_target\tcont_actual\tcontam_label\tcontam_source_acc\t"
               "n_contigs\ttotal_bp\n";
        for (const auto& r : rows)
            tsv << r.synth_acc << '\t' << r.source_acc << '\t' << r.rep << '\t'
                << r.comp_target << '\t' << r.comp_actual << '\t'
                << r.cont_target << '\t' << r.cont_actual << '\t'
                << r.contam_label << '\t' << r.contam_source_acc << '\t'
                << r.n_contigs << '\t' << r.total_bp_out << '\n';
    }
    {
        std::ofstream man(add_path);
        if (!man) throw std::runtime_error("sim: cannot write " + add_path.string());
        man << "accession\ttaxonomy\tfile_path\n";
        for (const auto& r : rows)
            man << r.synth_acc << '\t' << r.taxonomy << '\t'
                << fs::absolute(r.fa_path).string() << '\n';
    }

    spdlog::info("sim: wrote {} FASTA files → {} + {}",
                 rows.size(), gt_path.string(), add_path.string());
    return 0;
}

} // namespace genopack
