#include <bits/stdc++.h>
#include <gtl/phmap.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <sux/bits/EliasFano.hpp>

using namespace std;

// popcnt helper for 4 words (SIMD-style unrolled)
uint64_t popcnt4(uint64_t x) { return __builtin_popcountll(x); }
// SPIDER-inspired Elias-Fano: Fastest practical EF implementation
// - Interleaved high/low sampling for cache efficiency
// - SIMD popcnt for high-rank 
// - Binary search + sampling for low-rank (O(1))
// - From SEA 2025 "SPIDER: Improved Succinct Rank and Select" [web:40]
struct FastEliasFano {
    uint64_t n, U, l, lowMask;
    vector<uint64_t> low_packed;
    vector<uint64_t> high_packed;
    
    // Sampling (every 2^16 elements for O(1) rank)
    static constexpr uint64_t SAMPLE_RATE = 1ULL << 16;
    static constexpr uint64_t NSAMPLES = SAMPLE_RATE;
    vector<uint64_t> low_samples;  // low[i*NSAMPLES]
    vector<uint64_t> high_samples; // high-rank[i*NSAMPLES]
    
    FastEliasFano(uint64_t n_, uint64_t U_) : n(n_), U(U_) {
        long double ratio = (long double)U / max<uint64_t>(1, n);
        l = ratio > 1.0L ? min<uint32_t>(60u, (uint32_t)floorl(log2l(ratio))) : 0u;
        lowMask = (l == 64 ? ~0ULL : ((1ULL << l) - 1ULL));
    }
    
    void build_sorted(const vector<uint64_t>& positions) {
        assert(positions.size() == n);
        
        uint64_t low_bits = n * l;
        uint64_t high_bits = (U >> l) + n;
        
        low_packed.resize((low_bits + 63) / 64);
        high_packed.resize((high_bits + 63) / 64);
        low_samples.resize((n + NSAMPLES - 1) / NSAMPLES);
        high_samples.resize((n + NSAMPLES - 1) / NSAMPLES);
        
        // Build compressed structure
        uint64_t low_pos = 0, high_pos_max = 0;
        for (uint64_t i = 0; i < n; ++i) {
            uint64_t x = positions[i];
            uint64_t hi = x >> l;
            uint64_t lo = x & lowMask;
            
            // Pack low bits
            pack_low(low_pos, lo);
            low_pos += l;
            
            // Set high bit at hi + i
            uint64_t hpos = hi + i;
            high_packed[hpos >> 6] |= (1ULL << (hpos & 63));
            high_pos_max = max(high_pos_max, hpos + 1);
            
            // Sampling
            if (i % NSAMPLES == 0) {
                low_samples[i / NSAMPLES] = lo;
                high_samples[i / NSAMPLES] = rank_high_fast(hpos);
            }
        }
        
        // Finalize samples
        high_samples.back() = rank_high_fast(high_pos_max - 1);
    }
    
private:
    void pack_low(uint64_t bitpos, uint64_t value) {
        if (l == 0) return;
        uint64_t word = bitpos >> 6;
        uint32_t off = bitpos & 63;
        uint64_t v = value << off;
        low_packed[word] |= v;
        if (off + l > 64) {
            low_packed[word + 1] |= value >> (64 - off);
        }
    }
    
public:
    // SIMD-accelerated high-rank (SPIDER-style)
    uint64_t rank_high_fast(uint64_t pos) {
        if (pos >= (U >> l) + n - 1) return n;
        
        uint64_t word = pos >> 6;
        uint64_t bit = pos & 63;
        
        // Sample-accelerated: jump to nearby sample
        uint64_t sample_idx = word / (NSAMPLES * 64 / NSAMPLES);
        uint64_t base_rank = high_samples[sample_idx];
        uint64_t base_word = sample_idx * (NSAMPLES * 64 / NSAMPLES);
        
        // SIMD popcnt on 4-8 words at a time
        uint64_t cnt = base_rank;
        for (uint64_t w = base_word; w < word; w += 4) {
            if (w + 4 <= word) {
                cnt += popcnt4(high_packed[w + 0]);
                cnt += popcnt4(high_packed[w + 1]);
                cnt += popcnt4(high_packed[w + 2]);
                cnt += popcnt4(high_packed[w + 3]);
                w += 3; // adjust loop
            } else {
                cnt += __builtin_popcountll(high_packed[w]);
            }
        }
        
        // Final word partial
        cnt += __builtin_popcountll(high_packed[word] & ((1ULL << bit) - 1));
        return cnt;
    }
    
    // Extract i-th low value (SIMD extract)
    uint64_t get_low(uint64_t i) {
        uint64_t bitpos = i * l;
        uint64_t word = bitpos >> 6;
        uint32_t off = bitpos & 63;
        uint64_t v = (low_packed[word] >> off) & lowMask;
        if (off + l > 64) {
            v |= (low_packed[word + 1] << (64 - off)) & lowMask;
        }
        return v;
    }
    
    // Sampled low-rank: O(log n) binary search to O(1) with samples
    uint64_t rank_low_fast(uint64_t lo_target) {
        // Binary search between samples
        uint64_t left = 0, right = (n + NSAMPLES - 1) / NSAMPLES;
        while (left < right) {
            uint64_t mid = (left + right) / 2;
            if (low_samples[mid] <= lo_target) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        // Scan within sample block (max NSAMPLES=64K → fast)
        uint64_t base = left * NSAMPLES;
        uint64_t cnt = base;
        for (uint64_t i = base; i < min(base + NSAMPLES, n); ++i) {
            if (get_low(i) <= lo_target) cnt++;
            else break;
        }
        return cnt;
    }
    
    // Full O(1) rank1
    uint64_t rank1(uint64_t pos) {
        uint64_t hi = pos >> l;
        uint64_t lo = pos & lowMask;
        return rank_high_fast(hi) + rank_low_fast(lo);
    }
    
    double bytes() const {
        return (low_packed.size() + high_packed.size() + 
                low_samples.size() + high_samples.size()) * 8.0;
    }
};

void bench_ht() {

}


int main() {
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(1);

    const uint64_t N = 10'000'000ULL;
    const uint64_t U = 1ULL << (31+31);

    cout << "FastEliasFano (SPIDER): 10M 1s in [0,2^62)...\n";

    // 1. Generate sorted positions
    vector<uint64_t> positions;
    positions.reserve(N);
    
    mt19937_64 rng(123456789);
    long double expected_gap = (long double)U / N;
    uint64_t current = 0;
    
    auto build_start = chrono::high_resolution_clock::now();
    for (uint64_t i = 0; i < N; ++i) {
        long double u = (long double)rng() / ULLONG_MAX;
        if (u <= 0.0L) u = 1e-18L;
        long double gap = -logl(u) * expected_gap;
        uint64_t dg = (uint64_t)gap ?: 1;
        current += dg;
        if (current >= U || (i > 0 && current <= positions.back())) {
            current = positions.empty() ? 0 : positions.back() + 1;
        }
        positions.push_back(current);
    }

    cout << "no random kmers: " << positions.size() << "\n";
    
    // 2. Build Fast EF
    FastEliasFano ef(N, U);
    ef.build_sorted(positions);
    auto build_end = chrono::high_resolution_clock::now();
    
    cout << "✓ Built: " << ef.bytes()/1e6 << "MB (" << (ef.bytes()*8/N) << " bpe)\n";
    cout << "  Build time: " << chrono::duration<double>(build_end - build_start).count() << "s\n";

    // // 3. 100M random queries
    const uint64_t Q = 100'000'000ULL;
    vector<uint64_t> queries(Q);
    uniform_int_distribution<uint64_t> dist(0, U-1);
    for (auto& q : queries) q = dist(rng);

    // Warmup
    for (int i = 0; i < 100000; ++i) ef.rank1(dist(rng));

    auto query_start = chrono::high_resolution_clock::now();
    uint64_t total_rank = 0;
    for (uint64_t q : queries) {
        total_rank += ef.rank1(q+1)-ef.rank1(q);
    }
    auto query_end = chrono::high_resolution_clock::now();

    double ns_query = chrono::duration<double, nano>(query_end - query_start).count() / (Q);
    double Mqps = 1e3 / ns_query;

    cout << "\n✓ 100M random rank1: " << ns_query << " ns/query (" << Mqps << " Mqps)\n";
    cout << "  Total 1s: " << total_rank << "\n";


    cout << "bench SDSL\n";
    build_start = chrono::high_resolution_clock::now();
    using namespace seqan3::contrib::sdsl;
    sd_vector<> r;
    rank_support_sd<> r_rank;
    sd_vector_builder builder(U, N);
    for(uint64_t pos : positions)
        builder.set(pos);
    r = sd_vector<>(builder);
    r_rank = rank_support_sd<>(&r);
    build_end = chrono::high_resolution_clock::now();
    size_t r_size = size_in_bytes(r);
    cout << "✓ Built: " << r_size/1e6 << "MB (" << (r_size*8/N) << " bpe)\n";
    cout << "  Build time: " << chrono::duration<double>(build_end - build_start).count() << "s\n";

    // Warmup
    for (int i = 0; i < 100000; ++i) r[dist(rng)];

    query_start = chrono::high_resolution_clock::now();
    total_rank = 0;
    for (uint64_t q : queries) {
        total_rank += r[q];
    }
    query_end = chrono::high_resolution_clock::now();

    ns_query = chrono::duration<double, nano>(query_end - query_start).count() / (Q);
    Mqps = 1e3 / ns_query;

    cout << "\n✓ 100M random R: " << ns_query << " ns/query (" << Mqps << " Mqps)\n";
    cout << "  Total 1s: " << total_rank << "\n";


    cout << "bench SUX\n";
    build_start = chrono::high_resolution_clock::now();
    sux::bits::EliasFano r2 = sux::bits::EliasFano(positions, U);
    build_end = chrono::high_resolution_clock::now();
    uint64_t r2_size = r2.bitCount()/8;
    cout << "✓ Built: " << r2_size/1e6 << "MB (" << (r2_size*8/N) << " bpe)\n";
    cout << "  Build time: " << chrono::duration<double>(build_end - build_start).count() << "s\n";

    // Warmup
    for (int i = 0; i < 100000; ++i) {
        uint64_t q = dist(rng);
        bool e = r2.rank(q+1)-r2.rank(q);
    }
    query_start = chrono::high_resolution_clock::now();
    total_rank = 0;
    for (uint64_t q : queries) {
        total_rank += r2.rank(q+1)-r2.rank(q);
    }
    query_end = chrono::high_resolution_clock::now();

    ns_query = chrono::duration<double, nano>(query_end - query_start).count() / (Q);
    Mqps = 1e3 / ns_query;

    cout << "\n✓ 100M random R: " << ns_query << " ns/query (" << Mqps << " Mqps)\n";
    cout << "  Total 1s: " << total_rank << "\n";


    cout << "bench HT\n";
    build_start = chrono::high_resolution_clock::now();
    gtl::flat_hash_set<uint64_t> hashmap;
    for(uint64_t pos : positions) {
        hashmap.insert(pos);
    }
    build_end = chrono::high_resolution_clock::now();
    size_t ht_size = 65*hashmap.bucket_count()/8;
    cout << "✓ Built: " << ht_size/1e6 << "MB (" << (ht_size*8/N) << " bpe)\n";
    cout << "  Build time: " << chrono::duration<double>(build_end - build_start).count() << "s\n";

    // Warmup
    for (int i = 0; i < 100000; ++i) hashmap.contains(dist(rng));

    query_start = chrono::high_resolution_clock::now();
    total_rank = 0;
    for (uint64_t q : queries) {
        total_rank += hashmap.contains(q);
    }
    query_end = chrono::high_resolution_clock::now();

    ns_query = chrono::duration<double, nano>(query_end - query_start).count() / (Q);
    Mqps = 1e3 / ns_query;

    cout << "\n✓ 100M random HT: " << ns_query << " ns/query (" << Mqps << " Mqps)\n";
    cout << "  Total 1s: " << total_rank << "\n";


    return 0;
}

