#include <immintrin.h>
// ...existing code...


inline void UnitigsDictionaryHash2::fill_buffer(std::vector<uint64_t> &buffer, const uint64_t mask, const size_t p, const size_t q)
{
    // Precompute min(next_endpoint, o+k+k) for the range [p, q)
    std::vector<size_t> precomputed_endpoints(q - p);
    size_t total_elements = 0;
    for (size_t i = 0; i < q - p; ++i) {
        size_t o = offsets[p + i];
        size_t next_endpoint = endpoints_select(endpoints_rank(o + 1) + 1);
        size_t e = o + span + k;
        if (e > next_endpoint)
            e = next_endpoint;
        precomputed_endpoints[i] = e;
        total_elements += e - o;
    }

    buffer.resize(total_elements);
    size_t buf_idx = 0;

    constexpr int SIMD_WIDTH = 8; // AVX-512: 8 x 64-bit lanes
    size_t i = 0;
    for (; i + SIMD_WIDTH <= q - p; i += SIMD_WIDTH) {
        __m512i hash = _mm512_setzero_si512();
        size_t o[SIMD_WIDTH];
        size_t e[SIMD_WIDTH];

        for (int lane = 0; lane < SIMD_WIDTH; ++lane) {
            o[lane] = offsets[p + i + lane];
            e[lane] = precomputed_endpoints[i + lane];
        }

        // Compute initial hashes for each lane
        for (uint64_t j = 0; j < k; ++j) {
            uint64_t new_rank[SIMD_WIDTH];
            for (int lane = 0; lane < SIMD_WIDTH; ++lane) {
                new_rank[lane] = seqan3::to_rank(text[o[lane] + j]);
            }
            __m512i ranks = _mm512_loadu_si512(new_rank);
            hash = _mm512_slli_epi64(hash, 2);
            hash = _mm512_or_si512(hash, ranks);
            __m512i mask_vec = _mm512_set1_epi64(mask);
            hash = _mm512_and_si512(hash, mask_vec);
        }

        // Store initial hashes
        uint64_t hash_out[SIMD_WIDTH];
        _mm512_storeu_si512(hash_out, hash);
        for (int lane = 0; lane < SIMD_WIDTH; ++lane) {
            buffer[buf_idx++] = hash_out[lane];
        }

        // Continue for j = o[x]+k to e[x]
        for (uint64_t j = k; j < k + k; ++j) {
            uint64_t new_rank[SIMD_WIDTH];
            for (int lane = 0; lane < SIMD_WIDTH; ++lane) {
                if (o[lane] + j < e[lane])
                    new_rank[lane] = seqan3::to_rank(text[o[lane] + j]);
                else
                    new_rank[lane] = 0; // or some neutral value
            }
            __m512i ranks = _mm512_loadu_si512(new_rank);
            hash = _mm512_slli_epi64(hash, 2);
            hash = _mm512_or_si512(hash, ranks);
            __m512i mask_vec = _mm512_set1_epi64(mask);
            hash = _mm512_and_si512(hash, mask_vec);

            _mm512_storeu_si512(hash_out, hash);
            for (int lane = 0; lane < SIMD_WIDTH; ++lane) {
                if (o[lane] + j < e[lane])
                    buffer[buf_idx++] = hash_out[lane];
            }
        }
    }

    // Handle remaining hashes (less than SIMD_WIDTH)
    for (; i < q - p; ++i) {
        uint64_t hash = 0;
        size_t o = offsets[p + i];
        size_t e = precomputed_endpoints[i];
        for (uint64_t j = o; j < o + k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
        }
        buffer[buf_idx++] = hash;
        for (size_t j = o + k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buffer[buf_idx++] = hash;
        }
    }
}