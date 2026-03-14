#include <genopack/util.hpp>
#include <zlib.h>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace genopack {

std::string decompress_gz(const std::filesystem::path& path) {
    gzFile gz = gzopen(path.string().c_str(), "rb");
    if (!gz) {
        // Try as plain FASTA
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) throw std::runtime_error("Cannot open: " + path.string());
        size_t sz = static_cast<size_t>(f.tellg());
        f.seekg(0);
        std::string out(sz, '\0');
        f.read(out.data(), sz);
        return out;
    }
    std::string out;
    char buf[65536];
    int n;
    while ((n = gzread(gz, buf, sizeof(buf))) > 0)
        out.append(buf, n);
    gzclose(gz);
    return out;
}

uint64_t genome_minhash(const std::string& fasta, int k) {
    uint64_t min_hash = UINT64_MAX;
    uint64_t kmer = 0;
    int kmer_len = 0;
    for (char c : fasta) {
        if (c == '>' || c == '\n' || c == '\r') { kmer_len = 0; kmer = 0; continue; }
        char u = static_cast<char>(c & ~0x20);
        uint8_t base;
        if      (u == 'A') base = 0;
        else if (u == 'C') base = 1;
        else if (u == 'G') base = 2;
        else if (u == 'T') base = 3;
        else { kmer_len = 0; kmer = 0; continue; }
        kmer = (kmer << 2) | base;
        if (++kmer_len >= k) {
            uint64_t v = kmer & ((1ULL << (2 * k)) - 1);
            // MurmurHash3 finalizer mix
            v ^= v >> 33;
            v *= 0xff51afd7ed558ccdULL;
            v ^= v >> 33;
            v *= 0xc4ceb9fe1a85ec53ULL;
            v ^= v >> 33;
            if (v < min_hash) min_hash = v;
        }
    }
    return min_hash;
}

FastaStats compute_fasta_stats(const std::string& fasta, int k) {
    FastaStats s;
    uint64_t gc = 0;
    uint64_t at = 0;
    for (char c : fasta) {
        if (c == '>') {
            ++s.n_contigs;
        } else if (c != '\n' && c != '\r') {
            ++s.genome_length;
            char u = static_cast<char>(c & ~0x20);
            if      (u == 'G' || u == 'C') ++gc;
            else if (u == 'A' || u == 'T') ++at;
        }
    }
    s.oph_fingerprint = genome_minhash(fasta, k);
    uint64_t total = gc + at;
    s.gc_pct_x100 = (total > 0)
        ? static_cast<uint16_t>(gc * 10000 / total)
        : 0;
    return s;
}

uint32_t days_since_epoch() {
    using namespace std::chrono;
    auto now = system_clock::now();
    auto tp  = floor<days>(now);
    auto ref = sys_days{year{2024}/January/1};
    return static_cast<uint32_t>((tp - ref).count());
}

std::vector<BuildRecord> parse_tsv_records(const std::filesystem::path& tsv_path) {
    std::ifstream f(tsv_path);
    if (!f) throw std::runtime_error("Cannot open TSV: " + tsv_path.string());

    std::string header_line;
    std::getline(f, header_line);

    std::vector<std::string> cols;
    {
        std::istringstream ss(header_line);
        std::string tok;
        while (std::getline(ss, tok, '\t')) cols.push_back(tok);
    }

    auto find_col = [&](std::initializer_list<const char*> names) -> int {
        for (const char* name : names)
            for (int i = 0; i < (int)cols.size(); ++i)
                if (cols[i] == name) return i;
        return -1;
    };

    int idx_acc  = find_col({"accession", "acc"});
    int idx_path = find_col({"file_path", "path", "fasta_path", "fasta", "file"});
    int idx_comp = find_col({"completeness"});
    int idx_cont = find_col({"contamination"});
    int idx_glen = find_col({"genome_length"});
    int idx_ctg  = find_col({"n_contigs"});

    if (idx_acc  < 0) throw std::runtime_error("TSV missing 'accession' column");
    if (idx_path < 0) throw std::runtime_error("TSV missing file path column"
        " (expected: file_path, path, fasta_path, fasta, or file)");

    std::set<int> known;
    if (idx_acc  >= 0) known.insert(idx_acc);
    if (idx_path >= 0) known.insert(idx_path);
    if (idx_comp >= 0) known.insert(idx_comp);
    if (idx_cont >= 0) known.insert(idx_cont);
    if (idx_glen >= 0) known.insert(idx_glen);
    if (idx_ctg  >= 0) known.insert(idx_ctg);

    std::vector<int>         extra_indices;
    std::vector<std::string> extra_names;
    for (int i = 0; i < (int)cols.size(); ++i) {
        if (!known.count(i)) {
            extra_indices.push_back(i);
            extra_names.push_back(cols[i]);
        }
    }

    std::vector<BuildRecord> records;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> fields;
        {
            std::istringstream ss(line);
            std::string tok;
            while (std::getline(ss, tok, '\t')) fields.push_back(tok);
        }
        if ((int)fields.size() <= std::max(idx_acc, idx_path)) continue;

        BuildRecord r;
        r.accession = fields[idx_acc];
        r.file_path = fields[idx_path];
        if (idx_comp >= 0 && idx_comp < (int)fields.size())
            r.completeness  = std::stof(fields[idx_comp]);
        if (idx_cont >= 0 && idx_cont < (int)fields.size())
            r.contamination = std::stof(fields[idx_cont]);
        if (idx_glen >= 0 && idx_glen < (int)fields.size())
            r.genome_length  = std::stoull(fields[idx_glen]);
        if (idx_ctg  >= 0 && idx_ctg  < (int)fields.size())
            r.n_contigs      = static_cast<uint32_t>(std::stoul(fields[idx_ctg]));
        for (int j = 0; j < (int)extra_indices.size(); ++j) {
            int ci = extra_indices[j];
            r.extra_fields.emplace_back(extra_names[j],
                ci < (int)fields.size() ? fields[ci] : "");
        }
        records.push_back(std::move(r));
    }
    return records;
}

} // namespace genopack
