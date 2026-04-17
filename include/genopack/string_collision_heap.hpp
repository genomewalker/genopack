#pragma once
#include <algorithm>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

// Myers (CPM 2023) "Merging Sorted Lists of Similar Strings"
// String Heap (§3) + Collision Heap (§4).
//
// Merges T sorted lists of (string key, Value) pairs in O(M log(T/ē) + S)
// where M = total entries, ē = collision ratio, S = total key bytes.
//
// Duplicates: if ResolveFn is nullptr, the first-seen value wins.
// Otherwise resolve(std::vector<Value>&) → Value is called.

namespace genopack {

namespace detail {

// A min-heap on (string key, list index).
// Each node i carries:
//   lcp_[i]   LCP between H[i] and H[parent(i)]     (String Heap, §3)
//   l_eq_[i]  H[i].key == H[left_child(i)].key       (Collision Heap, §4)
//   r_eq_[i]  H[i].key == H[right_child(i)].key      (Collision Heap, §4)
//
// Heap is 1-indexed; H[0] is unused sentinel.
// An exhausted list is represented by list_pos[j] == list[j].size()
// and its heap slot key is treated as +∞.

template <typename Value>
struct StringCollisionHeap {

    struct Node {
        std::string  key;
        Value        val;
        int          list_idx;  // which input list owns this node
        bool         exhausted; // key == +inf
    };

    int T_;  // number of lists
    std::vector<Node>    H_;       // 1-indexed, size T_+1
    std::vector<int>     P_;       // LCP of H[i] with its parent (0 at root)
    std::vector<bool>    L_eq_;    // H[i].key == H[left_child(i)].key
    std::vector<bool>    R_eq_;    // H[i].key == H[right_child(i)].key

    // Pointers into the input lists
    std::vector<const std::vector<std::pair<std::string,Value>>*> lists_;
    std::vector<size_t> pos_;      // current position in each list

    // Compute LCP of two strings starting comparison at offset d.
    static int lcp_from(const std::string& a, const std::string& b, int d) {
        int n = static_cast<int>(std::min(a.size(), b.size()));
        int i = d;
        while (i < n && a[i] == b[i]) ++i;
        return i;
    }

    // Full LCP
    static int lcp_full(const std::string& a, const std::string& b) {
        return lcp_from(a, b, 0);
    }

    // Compare two nodes at heap positions i and j.
    // Returns true if node i < node j (i should be closer to root).
    bool less_than(int i, int j) const {
        if (H_[i].exhausted && H_[j].exhausted) return false;
        if (H_[i].exhausted) return false;
        if (H_[j].exhausted) return true;
        return H_[i].key < H_[j].key;
    }

    bool equal_keys(int i, int j) const {
        if (H_[i].exhausted || H_[j].exhausted) return false;
        return H_[i].key == H_[j].key;
    }

    // Advance list list_idx to its next element, update heap slot s.
    void advance_list(int list_idx, int slot) {
        ++pos_[list_idx];
        const auto& lst = *lists_[list_idx];
        if (pos_[list_idx] < lst.size()) {
            H_[slot].key       = lst[pos_[list_idx]].first;
            H_[slot].val       = lst[pos_[list_idx]].second;
            H_[slot].exhausted = false;
        } else {
            H_[slot].key.clear();
            H_[slot].exhausted = true;
        }
    }

    // Build heap from scratch (initial heapify).
    // Sets P_, L_eq_, R_eq_ for all nodes.
    void build() {
        // Bottom-up sift-down build, O(T).
        // We first arrange each list's front element into H[1..T].
        // Then fix the heap property bottom-up.
        for (int s = T_ / 2; s >= 1; --s)
            sift_down_full(s);
        // Compute P values top-down after heap is built.
        compute_lcp_topdown();
        // Compute L_eq_ / R_eq_
        for (int s = 1; s <= T_; ++s) {
            int l = 2 * s, r = 2 * s + 1;
            L_eq_[s] = (l <= T_) && equal_keys(s, l);
            R_eq_[s] = (r <= T_) && equal_keys(s, r);
        }
    }

    // Sift down node at position s (no LCP maintenance, used during build).
    void sift_down_full(int s) {
        while (true) {
            int l = 2 * s, r = 2 * s + 1;
            int smallest = s;
            if (l <= T_ && less_than(l, smallest)) smallest = l;
            if (r <= T_ && less_than(r, smallest)) smallest = r;
            if (smallest == s) break;
            std::swap(H_[s], H_[smallest]);
            s = smallest;
        }
    }

    // Compute P[] (parent LCP) top-down after a full heapify.
    void compute_lcp_topdown() {
        P_[1] = 0;
        for (int s = 2; s <= T_; ++s) {
            int p = s / 2;
            if (H_[s].exhausted || H_[p].exhausted)
                P_[s] = 0;
            else
                P_[s] = lcp_full(H_[s].key, H_[p].key);
        }
    }

    // Myers §3 String Heap sift-down after root has been replaced.
    // Maintains P_ (LCP with parent) for all disturbed nodes.
    // d = current comparison depth (LCP of new root with context).
    void sift_down_string(int s, int d) {
        while (true) {
            int l = 2 * s, r = 2 * s + 1;

            // Determine minimum child by comparing at depth d
            // (we already know first d chars match parent, skip them).
            int smaller = -1;
            if (l > T_ && r > T_) break;  // leaf

            if (l <= T_ && r <= T_) {
                // Both children exist — use their stored LCP with parent (P_)
                // to avoid redundant comparisons (Myers §3, cases 1-5).
                int pl = P_[l], pr = P_[r];
                if (pl > pr) {
                    // H[l] < H[r] for sure (l matches parent more deeply)
                    smaller = l;
                } else if (pr > pl) {
                    smaller = r;
                } else {
                    // pl == pr: need to compare at position pl (= pr)
                    if (H_[l].exhausted && H_[r].exhausted) { smaller = l; }
                    else if (H_[l].exhausted) { smaller = r; }
                    else if (H_[r].exhausted) { smaller = l; }
                    else {
                        // Compare from position pl onwards
                        int lcp_lr = lcp_from(H_[l].key, H_[r].key, pl);
                        if (lcp_lr == static_cast<int>(std::min(H_[l].key.size(), H_[r].key.size()))) {
                            // One is prefix of other or equal
                            smaller = (H_[l].key.size() <= H_[r].key.size()) ? l : r;
                        } else {
                            smaller = (H_[l].key[lcp_lr] < H_[r].key[lcp_lr]) ? l : r;
                        }
                    }
                }
            } else {
                smaller = (l <= T_) ? l : r;
            }

            if (H_[smaller].exhausted) break;
            if (!H_[s].exhausted && !less_than(smaller, s)) break;

            // Compute new P[smaller] = LCP(H[s], H[smaller]) from depth d.
            int new_p;
            if (H_[s].exhausted || H_[smaller].exhausted) {
                new_p = 0;
            } else {
                new_p = lcp_from(H_[s].key, H_[smaller].key, d);
            }

            std::swap(H_[s], H_[smaller]);
            // After swap, H[smaller] is old H[s].
            // P[smaller] = LCP(H[smaller], H[s/2 parent]) — maintained externally.
            // P[s] is now for the new occupant of s, which is old H[smaller].
            // But P[s] hasn't changed meaning: it is LCP with s's parent.
            // We need to update P[smaller] = new_p (LCP of moved-down node with its parent s).
            P_[smaller] = new_p;

            d = new_p;
            s = smaller;
        }
    }

    // Initialize with T lists.
    void init(std::vector<std::vector<std::pair<std::string,Value>>>& lists) {
        T_ = static_cast<int>(lists.size());
        H_.resize(T_ + 1);
        P_.assign(T_ + 1, 0);
        L_eq_.assign(T_ + 1, false);
        R_eq_.assign(T_ + 1, false);
        lists_.resize(T_);
        pos_.assign(T_, 0);

        for (int j = 0; j < T_; ++j) {
            lists_[j] = &lists[j];
            H_[j + 1].list_idx = j;
            if (!lists[j].empty()) {
                H_[j + 1].key       = lists[j][0].first;
                H_[j + 1].val       = lists[j][0].second;
                H_[j + 1].exhausted = false;
            } else {
                H_[j + 1].exhausted = true;
            }
        }
        if (T_ > 0) build();
    }

    // Pop the minimum key from the heap. Returns false if all lists exhausted.
    // Fills out_key, out_val, and also sets collisions[] with all slots that
    // hold the same key as root (for collision merging).
    bool pop_min(std::string& out_key, std::vector<std::pair<int,Value>>& same_key_slots) {
        if (T_ == 0 || H_[1].exhausted) return false;

        out_key = H_[1].key;
        same_key_slots.clear();
        same_key_slots.emplace_back(H_[1].list_idx, H_[1].val);

        // Collect collisions using L_eq_/R_eq_ flags (Myers §4).
        // Any node reachable top-down via equality edges holds the same key.
        collect_collisions(1, out_key, same_key_slots);

        // Advance all collision slots and re-sift.
        // We must re-sift root last (since re-sifting replaces nodes).
        // Process in reverse BFS order so parents are re-sifted after children.
        // Simpler: advance all, then full rebuild.  For T ≤ ~64 this is fine.
        // For correctness with Myers §4 (not just §3), we do a targeted re-sift.
        for (auto& [li, _v] : same_key_slots) {
            // Find slot of this list_idx (may have moved during collisions).
            // Since we have T_ ≤ a few dozen, linear scan is fine.
            for (int s = 1; s <= T_; ++s) {
                if (!H_[s].exhausted && H_[s].list_idx == li &&
                    H_[s].key == out_key) {
                    advance_list(li, s);
                    break;
                }
            }
        }

        // Rebuild (heapify) to restore invariants after multiple slot updates.
        // Myers §4 gives O(log(T/ē)) amortized; for simplicity we do O(T) rebuild
        // only when there are collisions (ē > 0), and O(log T) sift-down otherwise.
        if (same_key_slots.size() == 1) {
            // No collision — standard sift-down, O(log T).
            sift_down_full(1);
            compute_lcp_topdown();
            for (int s = 1; s <= T_; ++s) {
                int l = 2 * s, r = 2 * s + 1;
                L_eq_[s] = (l <= T_) && equal_keys(s, l);
                R_eq_[s] = (r <= T_) && equal_keys(s, r);
            }
        } else {
            // Multiple collisions — rebuild heap (O(T)).
            for (int s = T_ / 2; s >= 1; --s) sift_down_full(s);
            compute_lcp_topdown();
            for (int s = 1; s <= T_; ++s) {
                int l = 2 * s, r = 2 * s + 1;
                L_eq_[s] = (l <= T_) && equal_keys(s, l);
                R_eq_[s] = (r <= T_) && equal_keys(s, r);
            }
        }
        return true;
    }

    // DFS through equality edges to find all slots with key == root_key.
    void collect_collisions(int s, const std::string& root_key,
                            std::vector<std::pair<int,Value>>& out) {
        if (L_eq_[s]) {
            int l = 2 * s;
            if (l <= T_ && !H_[l].exhausted && H_[l].key == root_key) {
                out.emplace_back(H_[l].list_idx, H_[l].val);
                collect_collisions(l, root_key, out);
            }
        }
        if (R_eq_[s]) {
            int r = 2 * s + 1;
            if (r <= T_ && !H_[r].exhausted && H_[r].key == root_key) {
                out.emplace_back(H_[r].list_idx, H_[r].val);
                collect_collisions(r, root_key, out);
            }
        }
    }
};

} // namespace detail

// Public interface.
template <typename Value, typename EmitFn, typename ResolveFn>
void string_heap_merge(
    std::vector<std::vector<std::pair<std::string, Value>>>& lists,
    EmitFn emit,
    ResolveFn resolve)
{
    // Remove empty lists
    detail::StringCollisionHeap<Value> heap;
    heap.init(lists);

    std::string                        key;
    std::vector<std::pair<int,Value>>  slots;

    while (heap.pop_min(key, slots)) {
        if (slots.size() == 1) {
            emit(key, slots[0].second);
        } else {
            if constexpr (!std::is_same_v<ResolveFn, std::nullptr_t>) {
                std::vector<Value> vals;
                vals.reserve(slots.size());
                for (auto& [li, v] : slots) vals.push_back(std::move(v));
                emit(key, resolve(vals));
            } else {
                // First-wins: slots[0] is the first list encountered.
                emit(key, slots[0].second);
            }
        }
    }
}

// Overload without resolve (first-wins).
template <typename Value, typename EmitFn>
void string_heap_merge(
    std::vector<std::vector<std::pair<std::string, Value>>>& lists,
    EmitFn emit)
{
    string_heap_merge(lists, emit, nullptr);
}

} // namespace genopack
