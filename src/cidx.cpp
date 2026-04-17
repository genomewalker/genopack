#include <genopack/cidx.hpp>
#include <algorithm>
#include <stdexcept>
#include <thread>

namespace genopack {

// ── Hash ──────────────────────────────────────────────────────────────────────

uint64_t cidx_hash(std::string_view acc) {
    uint64_t h = 14695981039346656037ULL;
    for (unsigned char c : acc) {
        h ^= c;
        h *= 1099511628211ULL;
    }
    return h;
}

// ── FASTA parser ──────────────────────────────────────────────────────────────

void parse_fasta_contig_accessions(std::string_view fasta,
                                   const std::function<void(std::string_view)>& cb)
{
    size_t i = 0;
    const size_t n = fasta.size();
    while (i < n) {
        if (fasta[i] != '>') {
            while (i < n && fasta[i] != '\n') ++i;
            if (i < n) ++i;
            continue;
        }
        ++i; // skip '>'
        while (i < n && (fasta[i] == ' ' || fasta[i] == '\t')) ++i;
        size_t tok_start = i;
        while (i < n && fasta[i] != ' ' && fasta[i] != '\t' &&
               fasta[i] != '\n' && fasta[i] != '\r') ++i;
        if (i > tok_start)
            cb(fasta.substr(tok_start, i - tok_start));
        while (i < n && fasta[i] != '\n') ++i;
        if (i < n) ++i;
    }
}

// ── CidxWriter ────────────────────────────────────────────────────────────────

void CidxWriter::add(std::string_view contig_acc, uint32_t genome_id) {
    entries_.push_back({cidx_hash(contig_acc), genome_id, 0});
}

void CidxWriter::add_hash(uint64_t acc_hash, uint32_t genome_id) {
    entries_.push_back({acc_hash, genome_id, 0});
}

SectionDesc CidxWriter::finalize(AppendWriter& writer, uint64_t section_id, uint64_t batch_id) {
    std::sort(entries_.begin(), entries_.end(),
              [](const CidxEntry& a, const CidxEntry& b) {
                  return a.acc_hash < b.acc_hash;
              });

    CidxHeader hdr{};
    hdr.magic     = SEC_CIDX;
    hdr.version   = 1;
    hdr.n_entries = entries_.size();
    hdr.batch_id  = batch_id;

    uint64_t start = writer.current_offset();
    writer.append(&hdr, sizeof(hdr));
    if (!entries_.empty())
        writer.append(entries_.data(), entries_.size() * sizeof(CidxEntry));

    SectionDesc sd{};
    sd.type              = SEC_CIDX;
    sd.version           = 1;
    sd.section_id        = section_id;
    sd.file_offset       = start;
    sd.compressed_size   = sizeof(hdr) + entries_.size() * sizeof(CidxEntry);
    sd.uncompressed_size = sd.compressed_size;
    sd.item_count        = entries_.size();
    sd.aux0              = batch_id;
    return sd;
}

// ── CidxReader ────────────────────────────────────────────────────────────────

void CidxReader::open(const uint8_t* base, uint64_t offset, uint64_t size) {
    if (size < sizeof(CidxHeader))
        throw std::runtime_error("CIDX section too small");
    hdr_ = reinterpret_cast<const CidxHeader*>(base + offset);
    if (hdr_->magic != SEC_CIDX)
        throw std::runtime_error("CIDX bad magic");
    n_       = hdr_->n_entries;
    entries_ = reinterpret_cast<const CidxEntry*>(base + offset + sizeof(CidxHeader));
    if (sizeof(CidxHeader) + n_ * sizeof(CidxEntry) > size)
        throw std::runtime_error("CIDX section truncated");
}

uint32_t CidxReader::find(uint64_t acc_hash) const {
    if (n_ == 0) return UINT32_MAX;
    size_t lo = 0, hi = n_;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (entries_[mid].acc_hash < acc_hash) lo = mid + 1;
        else                                    hi = mid;
    }
    if (lo < n_ && entries_[lo].acc_hash == acc_hash)
        return entries_[lo].genome_id;
    return UINT32_MAX;
}

uint32_t CidxReader::find(std::string_view contig_acc) const {
    return find(cidx_hash(contig_acc));
}

uint64_t CidxReader::batch_id() const {
    return hdr_ ? hdr_->batch_id : 0;
}

void CidxReader::scan(const std::function<void(uint64_t, uint32_t)>& cb) const {
    for (uint64_t i = 0; i < n_; ++i)
        cb(entries_[i].acc_hash, entries_[i].genome_id);
}

// ── MergedCidxReader ──────────────────────────────────────────────────────────

void MergedCidxReader::add_section(const uint8_t* base, uint64_t offset, uint64_t size) {
    CidxReader r;
    r.open(base, offset, size);
    auto it = std::lower_bound(readers_.begin(), readers_.end(), r.batch_id(),
        [](const CidxReader& a, uint64_t bid) { return a.batch_id() > bid; });
    readers_.insert(it, std::move(r));
}

uint32_t MergedCidxReader::find_hash(uint64_t acc_hash) const {
    for (const auto& r : readers_) {
        uint32_t gid = r.find(acc_hash);
        if (gid != UINT32_MAX) return gid;
    }
    return UINT32_MAX;
}

uint32_t MergedCidxReader::find(std::string_view contig_acc) const {
    return find_hash(cidx_hash(contig_acc));
}

size_t MergedCidxReader::total_entries() const {
    size_t n = 0;
    for (const auto& r : readers_) n += r.n_entries();
    return n;
}

std::vector<CidxEntry> MergedCidxReader::merge_all() const {
    size_t total = 0;
    for (const auto& r : readers_) total += r.n_entries();

    std::vector<CidxEntry> all;
    all.reserve(total);
    for (const auto& r : readers_)
        r.scan([&](uint64_t h, uint32_t gid) { all.push_back({h, gid, 0}); });

    // stable_sort preserves the reader iteration order (newest-first from add_section)
    // so that dedup below always keeps the newest entry for each hash.
    std::stable_sort(all.begin(), all.end(),
                     [](const CidxEntry& a, const CidxEntry& b) { return a.acc_hash < b.acc_hash; });

    // Deduplicate — keep first occurrence (newest section iterated first)
    auto out = all.begin();
    for (auto it = all.begin(); it != all.end(); ++it) {
        if (out == all.begin() || it->acc_hash != (out - 1)->acc_hash)
            *out++ = *it;
    }
    all.erase(out, all.end());
    return all;
}

// ── Batch lookup ──────────────────────────────────────────────────────────────

// Prefetch distance in CidxEntry units (16 bytes each → 4 per cache line).
// 32 entries = 8 cache lines ahead — hides ~150ns DRAM latency at 5ns/entry pace.
static constexpr size_t CIDX_PF = 32;

void MergedCidxReader::merge_join_section(const CidxReader&      reader,
                                          const CidxQuery*       queries,
                                          size_t                 n_queries,
                                          uint32_t*              out_genome_ids,
                                          std::vector<CidxQuery>& unfound_out)
{
    const CidxEntry* entries   = reader.data();
    const size_t     n_entries = reader.n_entries();

    if (n_entries == 0 || n_queries == 0) {
        for (size_t i = 0; i < n_queries; ++i) unfound_out.push_back(queries[i]);
        return;
    }

    // Binary-search to find the entry range covered by [queries[0].hash, queries[n-1].hash].
    // Skips the bulk of the array when queries are clustered (e.g. one organism).
    const uint64_t q_min = queries[0].hash;
    const uint64_t q_max = queries[n_queries - 1].hash;

    // ei_start: first entry with acc_hash >= q_min
    size_t ei_start = 0;
    {
        size_t lo = 0, hi = n_entries;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (entries[mid].acc_hash < q_min) lo = mid + 1;
            else                                hi = mid;
        }
        ei_start = lo;
    }

    // ei_end: first entry with acc_hash > q_max (exclusive)
    size_t ei_end = n_entries;
    {
        size_t lo = ei_start, hi = n_entries;
        while (lo < hi) {
            size_t mid = lo + (hi - lo) / 2;
            if (entries[mid].acc_hash <= q_max) lo = mid + 1;
            else                                 hi = mid;
        }
        ei_end = lo;
    }

    // All queries fall outside this section's range — carry all forward
    if (ei_start == ei_end) {
        for (size_t i = 0; i < n_queries; ++i) unfound_out.push_back(queries[i]);
        return;
    }

    size_t qi = 0;        // query cursor
    size_t ei = ei_start; // entry cursor — starts at first relevant entry

    while (qi < n_queries && ei < ei_end) {
        // Prefetch entries and queries ahead
        if (ei + CIDX_PF < n_entries)
            __builtin_prefetch(entries + ei + CIDX_PF, 0, 1);
        if (qi + CIDX_PF < n_queries)
            __builtin_prefetch(queries + qi + CIDX_PF, 0, 1);

        const uint64_t qh = queries[qi].hash;
        const uint64_t eh = entries[ei].acc_hash;

        if (qh < eh) {
            // Query not in this section — carry forward
            unfound_out.push_back(queries[qi++]);
        } else if (qh > eh) {
            ++ei;
        } else {
            // Match — write result, advance query only
            // (multiple queries can share a hash; entry stays for next query)
            out_genome_ids[queries[qi].orig_idx] = entries[ei].genome_id;
            ++qi;
        }
    }
    // Drain remaining queries that fell past the end of this section
    while (qi < n_queries)
        unfound_out.push_back(queries[qi++]);
}

void MergedCidxReader::batch_find(const std::string_view* accs,
                                  uint32_t*               out_genome_ids,
                                  size_t                  n,
                                  size_t                  n_threads) const
{
    if (n == 0) return;
    std::fill(out_genome_ids, out_genome_ids + n, UINT32_MAX);
    if (readers_.empty()) return;

    if (n_threads == 0)
        n_threads = std::thread::hardware_concurrency();
    n_threads = std::max<size_t>(1, n_threads);

    // ── Phase 1: parallel hash ────────────────────────────────────────────────
    std::vector<CidxQuery> queries(n);

    if (n_threads == 1 || n < 4096) {
        for (size_t i = 0; i < n; ++i)
            queries[i] = {cidx_hash(accs[i]), static_cast<uint32_t>(i), 0};
    } else {
        const size_t chunk = (n + n_threads - 1) / n_threads;
        std::vector<std::thread> threads;
        threads.reserve(n_threads);
        for (size_t t = 0; t < n_threads; ++t) {
            const size_t start = t * chunk;
            const size_t end   = std::min(start + chunk, n);
            if (start >= end) break;
            threads.emplace_back([&, start, end] {
                for (size_t i = start; i < end; ++i)
                    queries[i] = {cidx_hash(accs[i]), static_cast<uint32_t>(i), 0};
            });
        }
        for (auto& th : threads) th.join();
    }

    // ── Phase 2: sort queries by hash ─────────────────────────────────────────
    std::sort(queries.begin(), queries.end(),
              [](const CidxQuery& a, const CidxQuery& b) { return a.hash < b.hash; });

    // ── Phase 3: merge-join per section, newest→oldest ────────────────────────
    // After each pass only unfound queries enter the next pass.
    std::vector<CidxQuery> unfound;
    unfound.reserve(n);

    const CidxQuery* cur_queries = queries.data();
    size_t           cur_n       = n;

    for (const auto& reader : readers_) {
        if (cur_n == 0) break;
        unfound.clear();
        merge_join_section(reader, cur_queries, cur_n, out_genome_ids, unfound);
        // Swap so next pass works on unfound only
        std::swap(queries, unfound);
        cur_queries = queries.data();
        cur_n       = queries.size();
    }
}

} // namespace genopack
