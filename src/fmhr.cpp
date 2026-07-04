#include <genopack/fmhr.hpp>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <unistd.h>

namespace genopack {

FmhrWriter::FmhrWriter(FmhrWriter&& o) noexcept
    : entries_(std::move(o.entries_)),
      spill_path_(std::move(o.spill_path_)),
      spill_fp_(o.spill_fp_),
      spill_index_(std::move(o.spill_index_))
{
    o.spill_fp_ = nullptr;
    o.spill_path_.clear();
}

FmhrWriter& FmhrWriter::operator=(FmhrWriter&& o) noexcept {
    if (this == &o) return *this;
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
    if (!spill_path_.empty()) { std::remove(spill_path_.c_str()); spill_path_.clear(); }
    entries_     = std::move(o.entries_);
    spill_path_  = std::move(o.spill_path_);
    spill_fp_    = o.spill_fp_;
    spill_index_ = std::move(o.spill_index_);
    o.spill_fp_ = nullptr;
    o.spill_path_.clear();
    return *this;
}

FmhrWriter::FmhrWriter(std::string spill_dir) {
    if (spill_dir.empty()) return;
    spill_path_ = spill_dir + "/fmhr_spill_" + std::to_string(::getpid()) + ".bin";
    spill_fp_ = std::fopen(spill_path_.c_str(), "wb");
    if (!spill_fp_)
        throw std::runtime_error("FmhrWriter: cannot open spill " + spill_path_);
}

FmhrWriter::~FmhrWriter() {
    if (spill_fp_) { std::fclose(spill_fp_); spill_fp_ = nullptr; }
    if (!spill_path_.empty()) std::remove(spill_path_.c_str());
}

void FmhrWriter::add(uint64_t genus_hash, std::vector<uint64_t> hashes) {
    if (!spill_fp_) {
        entries_.push_back({genus_hash, std::move(hashes)});
        return;
    }
    const int64_t offset = static_cast<int64_t>(::ftello(spill_fp_));
    const uint32_t n = static_cast<uint32_t>(hashes.size());
    std::fwrite(&genus_hash, sizeof(genus_hash), 1, spill_fp_);
    std::fwrite(&n, sizeof(n), 1, spill_fp_);
    if (n > 0) std::fwrite(hashes.data(), sizeof(uint64_t), n, spill_fp_);
    spill_index_.push_back({genus_hash, offset, n});
}

SectionDesc FmhrWriter::finalize(AppendWriter& w, uint64_t section_id,
                                 uint32_t k, uint32_t c) {
    if (spill_fp_) {
        std::fclose(spill_fp_); spill_fp_ = nullptr;

        // Sort index by (genus_hash, file_offset) so we read each genus's chunks
        // in file order (mostly sequential since input was genus-sorted).
        // Peak memory = one genus's deduplicated hash set at a time.
        std::stable_sort(spill_index_.begin(), spill_index_.end(),
            [](const SpillIndex& a, const SpillIndex& b){
                return a.genus_hash < b.genus_hash ||
                       (a.genus_hash == b.genus_hash && a.file_offset < b.file_offset);
            });

        FILE* rf = std::fopen(spill_path_.c_str(), "rb");
        if (!rf) throw std::runtime_error("FmhrWriter: cannot reopen spill " + spill_path_);

        std::vector<uint64_t> buf;
        size_t i = 0;
        while (i < spill_index_.size()) {
            const uint64_t gh = spill_index_[i].genus_hash;
            std::vector<uint64_t> all_hashes;
            // Gather all chunks for this genus in file order.
            while (i < spill_index_.size() && spill_index_[i].genus_hash == gh) {
                const auto& si = spill_index_[i];
                // Record: [genus_hash:8][n_hashes:4][hashes:n*8]; hashes at +12.
                if (::fseeko(rf, static_cast<off_t>(si.file_offset) + 12, SEEK_SET) != 0)
                    throw std::runtime_error("FmhrWriter: seek failed");
                const size_t base = all_hashes.size();
                all_hashes.resize(base + si.n_hashes);
                if (si.n_hashes > 0 &&
                    std::fread(all_hashes.data() + base, sizeof(uint64_t), si.n_hashes, rf)
                        != si.n_hashes)
                    throw std::runtime_error("FmhrWriter: read failed");
                ++i;
            }
            if (!all_hashes.empty()) {
                std::sort(all_hashes.begin(), all_hashes.end());
                all_hashes.erase(std::unique(all_hashes.begin(), all_hashes.end()),
                                 all_hashes.end());
            }
            entries_.push_back({gh, std::move(all_hashes)});
        }
        std::fclose(rf);
        std::remove(spill_path_.c_str());
        spill_path_.clear();
        spill_index_.clear();
    }

    const uint32_t n = static_cast<uint32_t>(entries_.size());
    const uint32_t min_buckets = (n == 0) ? 1
        : static_cast<uint32_t>(static_cast<double>(n) / 0.7 + 1.0);
    const uint32_t n_buckets = next_pow2(min_buckets);
    const uint32_t mask      = n_buckets - 1;

    std::vector<uint32_t> buckets(n_buckets, FMHR_EMPTY);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t slot = static_cast<uint32_t>(entries_[i].genus_hash) & mask;
        while (buckets[slot] != FMHR_EMPTY) slot = (slot + 1) & mask;
        buckets[slot] = i;
    }

    w.align(8);
    const uint64_t section_start = w.current_offset();

    const uint64_t entries_off = sizeof(FmhrHeader);
    const uint64_t buckets_off = entries_off + static_cast<uint64_t>(n) * sizeof(FmhrEntry);
    const uint64_t hashes_off_unaligned = buckets_off + static_cast<uint64_t>(n_buckets) * 4;
    const uint64_t hashes_off = (hashes_off_unaligned + 7) & ~uint64_t{7};

    std::vector<FmhrEntry> disk_entries(n);
    uint64_t cur_hash_off = hashes_off;
    for (uint32_t i = 0; i < n; ++i) {
        disk_entries[i].genus_hash    = entries_[i].genus_hash;
        disk_entries[i].hashes_offset = cur_hash_off;
        disk_entries[i].n_hashes      = static_cast<uint32_t>(entries_[i].hashes.size());
        disk_entries[i]._pad          = 0;
        cur_hash_off += static_cast<uint64_t>(entries_[i].hashes.size()) * sizeof(uint64_t);
    }

    FmhrHeader hdr{};
    hdr.magic          = SEC_FMHR;
    hdr.n_genera       = n;
    hdr.n_buckets      = n_buckets;
    hdr.k              = k;
    hdr.c              = c;
    hdr.entries_offset = entries_off;
    hdr.buckets_offset = buckets_off;
    hdr.hashes_offset  = hashes_off;

    w.append(&hdr, sizeof(hdr));
    if (n > 0)
        w.append(disk_entries.data(), static_cast<uint64_t>(n) * sizeof(FmhrEntry));
    w.append(buckets.data(), static_cast<uint64_t>(n_buckets) * sizeof(uint32_t));

    if (const uint64_t pad = hashes_off - (buckets_off + n_buckets * 4); pad > 0) {
        uint64_t zero = 0;
        w.append(&zero, pad);
    }

    for (uint32_t i = 0; i < n; ++i) {
        if (!entries_[i].hashes.empty())
            w.append(entries_[i].hashes.data(),
                     static_cast<uint64_t>(entries_[i].hashes.size()) * sizeof(uint64_t));
    }

    const uint64_t section_end = w.current_offset();

    SectionDesc sd{};
    sd.type              = SEC_FMHR;
    sd.version           = 1;
    sd.flags             = 0;
    sd.section_id        = section_id;
    sd.file_offset       = section_start;
    sd.compressed_size   = section_end - section_start;
    sd.uncompressed_size = section_end - section_start;
    sd.item_count        = n;
    return sd;
}

} // namespace genopack
