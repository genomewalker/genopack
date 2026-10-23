// ArchiveReader::open must read the CATL, ACCX, GIDX, QCOL and TAXN sections with
// read syscalls (prefetch) instead of leaving them to MADV_RANDOM mmap faults. rchar
// counts read syscalls only, so it must grow by every byte of those sections. The
// fixture is built without GSTX/GCOV/FCOV, which open preads itself, so the other
// open reads stay smaller than any one section and a dropped section shows up.
#include <genopack/archive.hpp>
#include <genopack/format.hpp>
#include <genopack/mmap_file.hpp>
#include <genopack/toc.hpp>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include "test_utils.hpp"

using namespace genopack;

static uint64_t rchar() {
    std::ifstream io("/proc/self/io");
    std::string key;
    uint64_t value;
    while (io >> key >> value) if (key == "rchar:") return value;
    throw std::runtime_error("no rchar in /proc/self/io");
}

int main() {
    using namespace genopack_test;
    TempDir tmp("genopack_open_prefetch");
    {
        std::ofstream tsv(tmp.path / "in.tsv");
        tsv << "accession\ttaxonomy\tfile_path\n";
        for (int i = 0; i < 1500; ++i) {
            const auto fa = tmp.path / ("g" + std::to_string(i) + ".fa");
            write_fasta(fa, "g" + std::to_string(i), make_random_sequence(i + 1, 2000));
            tsv << "GCA_" << 100000000 + i << ".1\td__Bacteria;p__P" << i % 7 << ";c__C" << i % 13
                << ";o__O" << i % 29 << ";f__F" << i % 61 << ";g__G" << i % 151 << ";s__G" << i % 151
                << " sp" << i % 300 << "\t" << fa.string() << "\n";
        }
    }
    const auto gpk = tmp.path / "a.gpk";
    const std::string log = " >> " + shell_quote((tmp.path / "build.log").string()) + " 2>&1";
    run_checked(shell_quote(GENOPACK_BIN) + " build --no-gstx -i " + shell_quote((tmp.path / "in.tsv").string()) +
                " -o " + shell_quote(gpk.string()) + log);
    run_checked(shell_quote(GENOPACK_BIN) + " check " + shell_quote(gpk.string()) + log); // adds QCOL

    std::map<uint32_t, uint64_t> bytes;
    {
        MmapFileReader m;
        m.open(gpk);
        for (const auto& sd : TocReader::read(m).sections) bytes[sd.type] += sd.compressed_size;
    }
    uint64_t wanted = 0, smallest = UINT64_MAX;
    for (uint32_t type : {SEC_CATL, SEC_ACCX, SEC_GIDX, SEC_QCOL, SEC_TAXN}) {
        require(bytes[type] > 0, "fixture lacks a prefetched section type");
        wanted += bytes[type];
        smallest = std::min(smallest, bytes[type]);
    }

    const uint64_t before = rchar();
    ArchiveReader ar;
    ar.open(gpk);
    const uint64_t read = rchar() - before;
    std::cerr << "open read " << read << " bytes; prefetched sections " << wanted
              << " bytes, smallest " << smallest << "\n";
    require(read >= wanted, "open did not read every prefetched section byte");
    require(read - wanted < smallest, "fixture: other open reads could hide a dropped section");
    require(ar.genome_meta_by_accession("GCA_100000007.1").has_value(), "accession lookup failed after open");
    std::cout << "PASS\n";
}
