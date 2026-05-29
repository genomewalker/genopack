#pragma once
#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <string_view>

namespace genopack {

struct TnfTables {
    std::array<uint8_t, 256> base2bit;
    std::array<uint8_t, 256> canon;

    TnfTables() {
        base2bit.fill(0xFF);
        for (char c : {'A','a'}) base2bit[static_cast<uint8_t>(c)] = 0;
        for (char c : {'C','c'}) base2bit[static_cast<uint8_t>(c)] = 1;
        for (char c : {'G','g'}) base2bit[static_cast<uint8_t>(c)] = 2;
        for (char c : {'T','t'}) base2bit[static_cast<uint8_t>(c)] = 3;

        std::array<int, 256> canonical{};
        canonical.fill(-1);
        int next = 0;
        for (int i = 0; i < 256; ++i) {
            if (canonical[i] >= 0) continue;
            int b0 = (i >> 6) & 3, b1 = (i >> 4) & 3,
                b2 = (i >> 2) & 3, b3 = i & 3;
            int rc = (3-b3)*64 + (3-b2)*16 + (3-b1)*4 + (3-b0);
            canonical[i] = canonical[rc] = next++;
            canon[i] = canon[rc] = static_cast<uint8_t>(canonical[i]);
        }
    }
};

inline const TnfTables& tnf_tables() noexcept {
    static const TnfTables t;
    return t;
}

inline bool compute_tnf(std::string_view seq, Eigen::VectorXf& out) {
    if (seq.size() < 5000) return false;
    out = Eigen::VectorXf::Zero(136);

    const auto& tbl = tnf_tables();
    const uint8_t* b2b   = tbl.base2bit.data();
    const uint8_t* canon = tbl.canon.data();
    uint32_t counts[256] = {};

    uint32_t idx = 0, valid = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
        uint8_t bits = b2b[static_cast<uint8_t>(seq[i])];
        if (bits == 0xFF) { valid = 0; idx = 0; continue; }
        idx = ((idx << 2) | bits) & 0xFF;
        if (++valid >= 4) ++counts[idx];
    }

    for (int i = 0; i < 256; ++i)
        out[canon[i]] += static_cast<float>(counts[i]);

    float total = out.sum();
    if (total < 100.f) { out.setZero(); return false; }
    out /= total;
    float norm = out.norm();
    if (norm < 1e-8f) return false;
    out /= norm;
    return true;
}

// 16-dim Karlin ρ* (dinucleotide relative abundance).
// rho[i*4+j] = P(XiXj) / (P(Xi)*P(Xj)), ordering A=0 C=1 G=2 T=3.
// Returns false if sequence has fewer than 1000 valid bases.
inline bool compute_rho(std::string_view seq, float rho[16]) noexcept {
    uint32_t mono[4] = {}, di[16] = {};
    int prev = -1;
    for (char c : seq) {
        int b;
        switch (c | 0x20) {
            case 'a': b = 0; break; case 'c': b = 1; break;
            case 'g': b = 2; break; case 't': b = 3; break;
            default: prev = -1; continue;
        }
        mono[b]++;
        if (prev >= 0) di[prev*4+b]++;
        prev = b;
    }
    uint32_t n_mono = mono[0]+mono[1]+mono[2]+mono[3];
    uint32_t n_di = 0; for (int i=0;i<16;i++) n_di += di[i];
    if (n_mono < 1000 || n_di == 0) return false;
    float fm[4], fd[16];
    for (int i=0;i<4;i++) fm[i] = static_cast<float>(mono[i]) / n_mono;
    for (int i=0;i<16;i++) fd[i] = static_cast<float>(di[i]) / n_di;
    for (int i=0;i<4;i++)
        for (int j=0;j<4;j++) {
            float den = fm[i]*fm[j];
            rho[i*4+j] = (den > 1e-10f) ? fd[i*4+j] / den : 1.0f;
        }
    return true;
}

} // namespace genopack
