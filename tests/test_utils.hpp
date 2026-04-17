#pragma once
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace genopack_test {

struct TempDir {
    std::filesystem::path path;

    explicit TempDir(const std::string& tag) {
        path = std::filesystem::temp_directory_path() /
            (tag + "_" + std::to_string(::getpid()));
        std::filesystem::remove_all(path);
        std::filesystem::create_directories(path);
    }

    ~TempDir() {
        std::error_code ec;
        std::filesystem::remove_all(path, ec);
    }
};

inline std::string make_sequence(const std::string& motif, size_t n) {
    std::string out;
    out.reserve(n);
    while (out.size() < n) out += motif;
    out.resize(n);
    return out;
}

inline std::string mutate_every(const std::string& input, size_t step) {
    std::string out = input;
    for (size_t i = step; i < out.size(); i += step)
        out[i] = (out[i] == 'A') ? 'T' : 'A';
    return out;
}

inline void write_fasta(const std::filesystem::path& path,
                        const std::string& header,
                        const std::string& seq) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write fasta");
    out << ">" << header << "\n";
    for (size_t i = 0; i < seq.size(); i += 60)
        out << seq.substr(i, 60) << "\n";
}

inline void write_tsv(const std::filesystem::path& path,
                      const std::vector<std::pair<std::string, std::filesystem::path>>& rows) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write tsv");
    out << "accession\tfile_path\n";
    for (const auto& [acc, fasta] : rows)
        out << acc << "\t" << fasta.string() << "\n";
}

inline int run_checked(const std::string& cmd) {
    int rc = std::system(cmd.c_str());
    if (rc != 0)
        throw std::runtime_error("command failed: " + cmd);
    return rc;
}

inline std::string read_text(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot read file");
    return std::string((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
}

inline std::string shell_quote(const std::string& s) {
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') out += "'\\''";
        else out.push_back(c);
    }
    out.push_back('\'');
    return out;
}

inline void require(bool cond, const std::string& msg) {
    if (!cond) throw std::runtime_error(msg);
}

} // namespace genopack_test
