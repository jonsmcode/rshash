#include <filesystem>
#include <immintrin.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rsindex3_simd.hpp"
#include "io.hpp"
#include "minimiser_rev_xor_views3.hpp"


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
}


RSIndexComp::RSIndexComp() : endpoints(std::vector<uint64_t>{}, 1) {}

RSIndexComp::RSIndexComp(
    uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
    uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3)
    : k(k), m1(m1), m2(m2), m3(m3),
      m_thres1(m_thres1), m_thres2(m_thres2), m_thres3(m_thres3),
      endpoints(std::vector<uint64_t>{}, 1)
{}


int RSIndexComp::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});

    std::cout << +m1 << " " << +m2 << " " << +m3 << "\n";

    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);
    const uint64_t M3 = 1ULL << (m3+m3);

    std::cout << "find minimizers...\n";
    bit_vector r1tmp = bit_vector(M1, 0);

    size_t N = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view1) {
            r1tmp[minimiser.minimiser_value] = 1;
        }
        N += record.size();
        no_sequences++;
    }
    rank_support_v<1> r1tmp_rank = rank_support_v<1>(&r1tmp);

    std::cout << "get sequences...\n";
    sequences = bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), N+1);
    // endpoints_ = seqan3::contrib::sdsl::sd_vector<>(sequences_);
    // endpoints_rank = seqan3::contrib::sdsl::rank_support_sd<>(&endpoints);
    // endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);

    std::cout << "count minimizers...\n";
    size_t c1tmp = r1tmp_rank(M1);
    uint8_t* count1tmp = new uint8_t[c1tmp];
    std::memset(count1tmp, 0, c1tmp*sizeof(uint8_t));

    uint64_t kmers = 0;
    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t i = r1tmp_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;
            size_t w = o/span + 1;
            count1tmp[i] += w;
            if(count1tmp[i] > m_thres1)
                count1tmp[i] = m_thres1;
            kmers += o;
            n += w;
        }
    }

    std::cout << "filling R_1...\n";
    r1 = bit_vector(M1, 0);
    for(auto & sequence : input) {
        for(auto && minimisers : sequence | view1) {
            if(count1tmp[r1tmp_rank(minimisers.minimiser_value)] < m_thres1)
                r1[minimisers.minimiser_value] = 1;
        }
    }
    r1_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r1);

    delete[] count1tmp;

    std::cout << "count minimizers1 again...\n";
    size_t c1 = r1_rank(M1);
    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                count1[i] += minimiser.occurrences/span + 1;
            }
        }
    }
    uint64_t n1 = 0;
    for(size_t i=0; i < c1; i++)
        n1 += count1[i];

    std::cout << "filling bitvector S_1...\n";
    s1 = bit_vector(n1+1, 0);
    s1[0] = 1;
    j = 0;
    for (size_t i=0; i < c1; i++) {
        j += count1[i];
        s1[j] = 1;
    }
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), n1+1, 3);

    std::cout << "filling offsets_1...\n";
    const size_t offset_width = std::bit_width(N);
    offsets1.width(offset_width);
    offsets1.resize(n1);

    std::memset(count1, 0, c1*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                size_t s = s1_select.select(i);
                size_t o = minimiser.occurrences;
                uint64_t j = 0;
                while(o > span) {
                    offsets1[s + count1[i]] = length + minimiser.range_position + j*span;
                    count1[i]++;
                    o -= span;
                    j++;
                }
                offsets1[s + count1[i]] = length + minimiser.range_position + j*span;
                count1[i]++;
            }
        }
        length += sequence.size();
    }

    delete[] count1;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers1;
    std::vector<size_t> skmer_positions;
    length = 0;
    for(auto & sequence : input) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view1) {
            level_up = r1[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view1) {
            current_level_up = r1[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers1.push_back(skmer);
                skmer_positions.push_back(length + start_position);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers1.push_back(skmer);
            skmer_positions.push_back(length + start_position);
        }

        length += sequence.size();
    }
    
    size_t len_rem_seqs = 0;
    for(auto & sequence : freq_skmers1)
        len_rem_seqs += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers1.size() << " (" << (double) freq_skmers1.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";


    std::cout << "build level 2...\n";
    auto view2 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    bit_vector r2tmp = bit_vector(M2, 0);

    for(auto & record : freq_skmers1) {
        for(auto && minimiser : record | view2)
            r2tmp[minimiser.minimiser_value] = 1;
    }
    rank_support_v<1> r2tmp_rank = rank_support_v<1>(&r2tmp);

    std::cout << "count minimizers2...\n";
    size_t c2tmp = r2tmp_rank(M2);
    uint8_t* count2tmp = new uint8_t[c2tmp];
    std::memset(count2tmp, 0, c2tmp*sizeof(uint8_t));

    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            size_t i = r2tmp_rank(minimiser.minimiser_value);
            count2tmp[i] += minimiser.occurrences/span + 1;
            if(count2tmp[i] > m_thres2)
                count2tmp[i] = m_thres2;
        }
    }

    std::cout << "fill R2...\n";
    bit_vector r2_ = bit_vector(M2, 0);
    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            if(count2tmp[r2tmp_rank(minimiser.minimiser_value)] < m_thres2)
                r2_[minimiser.minimiser_value] = 1;
        }
    }
    r2 = sd_vector<>(r2_);
    r2_rank = rank_support_sd<>(&r2);

    std::cout << "fill count 2...\n";
    size_t c2 = r2_rank(M2);
    uint8_t* count2 = new uint8_t[c2];
    std::memset(count2, 0, c2*sizeof(uint8_t));

    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view2) {
            count2[r2_rank(minimiser.minimiser_value)] = count2tmp[r2tmp_rank(minimiser.minimiser_value)];
        }
    }

    delete[] count2tmp;

    std::cout << "filling bitvector S_2...\n";
    uint64_t n2 = 0;
    for(uint64_t i=0; i < c2; i++) {
        n2 += count2[i];
    }

    s2 = bit_vector(n2+1, 0);
    s2[0] = 1;
    j = 0;
    for (size_t i=0; i < c2; i++) {
        j += count2[i];
        s2[j] = 1;
    }
    s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), n2+1, 3);

    std::cout << "filling offsets_2...\n";
    offsets2.width(offset_width);
    offsets2.resize(n2);
    std::memset(count2, 0, c2*sizeof(uint8_t));

    uint64_t skmer_idx = 0;
    for(auto & skmer : freq_skmers1) {
        for (auto && minimiser : skmer | view2) {
            if(r2[minimiser.minimiser_value]) {
                size_t i = r2_rank(minimiser.minimiser_value);
                size_t s = s2_select.select(i);
                size_t o = minimiser.occurrences;
                uint64_t j = 0;
                while(o > span) {
                    offsets2[s + count2[i]] = skmer_positions[skmer_idx] + minimiser.range_position + j*span;
                    count2[i]++;
                    o -= span;
                    j++;
                }
                offsets2[s + count2[i]] = skmer_positions[skmer_idx] + minimiser.range_position + j*span;
                count2[i]++;
            }
        }
        skmer_idx++;
    }

    delete[] count2;


    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers2;
    std::vector<size_t> skmer_positions2;
    skmer_idx = 0; 
    for(auto & sequence : freq_skmers1) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view2) {
            level_up = r2[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view2) {
            current_level_up = r2[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers2.push_back(skmer);
                skmer_positions2.push_back(skmer_positions[skmer_idx] + start_position);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers2.push_back(skmer);
            skmer_positions2.push_back(skmer_positions[skmer_idx] + start_position);
        }

        skmer_idx++;
    }
    
    len_rem_seqs = 0;
    for(auto & sequence : freq_skmers2)
        len_rem_seqs += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers2.size() << " (" << (double) freq_skmers2.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";


    std::cout << "build level 3...\n";

    auto view3 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    bit_vector r3tmp = bit_vector(M3, 0);

    for(auto & record : freq_skmers2) {
        for(auto && minimiser : record | view3)
            r3tmp[minimiser.minimiser_value] = 1;
    }
    rank_support_v<1> r3tmp_rank = rank_support_v<1>(&r3tmp);

    std::cout << "count minimizers3...\n";
    size_t c3tmp = r3tmp_rank(M3);
    uint16_t* count3tmp = new uint16_t[c3tmp];
    std::memset(count3tmp, 0, c3tmp*sizeof(uint16_t));

    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            size_t i = r3tmp_rank(minimiser.minimiser_value);
            count3tmp[i] += minimiser.occurrences/span + 1;
            if(count3tmp[i] > m_thres3)
                count3tmp[i] = m_thres3;
        }
    }

    std::cout << "fill R3...\n";
    bit_vector r3_ = bit_vector(M3, 0);
    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            if(count3tmp[r3tmp_rank(minimiser.minimiser_value)] < m_thres3)
                r3_[minimiser.minimiser_value] = 1;
        }
    }
    r3 = sd_vector<>(r3_);
    r3_rank = rank_support_sd<>(&r3);

    std::cout << "fill count 3...\n";
    size_t c3 = r3_rank(M3);
    uint16_t* count3 = new uint16_t[c3];
    std::memset(count3, 0, c3*sizeof(uint16_t));

    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | view3) {
            count3[r3_rank(minimiser.minimiser_value)] = count3tmp[r3tmp_rank(minimiser.minimiser_value)];
        }
    }

    delete[] count3tmp;

    std::cout << "filling bitvector S_3...\n";
    uint64_t n3 = 0;
    for(size_t i=0; i < c3; i++) {
        n3 += count3[i];
    }

    s3 = bit_vector(n3+1, 0);
    s3[0] = 1;
    j = 0;
    for (size_t i=0; i < c3; i++) {
        j += count3[i];
        s3[j] = 1;
    }
    s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), n3+1, 3);

    std::cout << "filling offsets_3...\n";
    offsets3.width(offset_width);
    offsets3.resize(n3);
    std::memset(count3, 0, c3*sizeof(uint16_t));

    skmer_idx = 0;
    for(auto & skmer : freq_skmers2) {
        for (auto && minimiser : skmer | view3) {
            if(r3[minimiser.minimiser_value]) {
                size_t i = r3_rank(minimiser.minimiser_value);
                size_t s = s3_select.select(i);
                size_t o = minimiser.occurrences;
                uint64_t j = 0;
                while(o > span) {
                    offsets3[s + count3[i]] = skmer_positions2[skmer_idx] + minimiser.range_position + j*span;
                    count3[i]++;
                    o -= span;
                    j++;
                }
                offsets3[s + count3[i]] = skmer_positions2[skmer_idx] + minimiser.range_position + j*span;
                count3[i]++;
            }
        }
        skmer_idx++;
    }

    delete[] count3;

    std::cout << "build level 4, HT...\n";

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers3;
    for(auto & sequence : freq_skmers2) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view3) {
            level_up = r3[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view3) {
            current_level_up = r3[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers3.push_back(skmer);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers3.push_back(skmer);
        }
    }

    len_rem_seqs = 0;
    for(auto & sequence : freq_skmers3)
        len_rem_seqs += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers3.size() << " (" << (double) freq_skmers3.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";

    std::cout << "filling HT...\n";
    auto view4 = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});
    for(auto & sequence : freq_skmers3) {
        for(auto && minimiser : sequence | view4) {
            hashmap.insert(std::min<uint64_t>(minimiser.window_value, minimiser.window_value_rev));
        }
    }

    std::cout << "copy text...\n";
    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << c1tmp << "\n";
    std::cout << "minimiser going to level 2: " << c1tmp-c1 << "  " << (double) (c1tmp-c1)/c1tmp*100 << "%\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1 << "\n";
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    std::cout << "no minimiser2: " << n2 << "\n";
    std::cout << "no distinct minimiser2: " << c2 << "\n";
    std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';
    std::cout << "minimiser going to level 3: " << c2tmp-c2 << "  " << (double) (c2tmp-c2)/c2tmp*100 << "%\n";
    std::cout << "no minimiser3: " << n3 << "\n";
    std::cout << "no distinct minimiser3: " << c3 << "\n";
    std::cout << "avg superkmers3: " << (double) n3/c3 <<  '\n';
    std::cout << "minimiser going to level 4: " << c3tmp-c3 << "  " << (double) (c3tmp-c3)/c3tmp*100 << "%\n";
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/kmers*100 << "%\n";

    std::cout << "density r1: " << (double) c1/M1*100 << "%\n";
    std::cout << "density r2: " << (double) c2/M2*100 << "%\n";
    std::cout << "density r3: " << (double) c3/M3*100 << "%\n";
    std::cout << "density s1: " << (double) s1_select.bitCount()/(n1+1)*100 <<  "%\n";
    std::cout << "density s2: " << (double) s2_select.bitCount()/(n2+1)*100 <<  "%\n";
    std::cout << "density s3: " << (double) s3_select.bitCount()/(n3+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/kmers << "\n";
    std::cout << "offsets1: " << (double) n1*offset_width/kmers << "\n";
    std::cout << "offsets2: " << (double) n2*offset_width/kmers << "\n";
    std::cout << "offsets3: " << (double) n3*offset_width/kmers << "\n";
    std::cout << "Hashtable: " << (double) 64*hashmap.size()/kmers << "\n";
    std::cout << "R_1: " << (double) M1/kmers << "\n";
    std::cout << "R_2: " << (double) 8*size_in_bytes(r2)/kmers << "\n";
    std::cout << "R_3: " << (double) 8*size_in_bytes(r3)/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1)/kmers << "\n";
    std::cout << "S_2: " << (double) (n2+1)/kmers << "\n";
    std::cout << "S_3: " << (double) (n3+1)/kmers << "\n";

    std::cout << "total: " << (double) (n1*offset_width+n2*offset_width+n3*offset_width+2*N+M1+8*size_in_bytes(r2)+8*size_in_bytes(r3)+n1+1+n2+1+n3+1+endpoints.bitCount()+64*hashmap.size())/kmers << "\n";

    return 0;
}


template<int level>
inline void RSIndexComp::fill_buffer_avx512(std::vector<uint64_t> &buffer, const uint64_t mask, size_t p, size_t q)
{
    size_t i = 0;
    for (; i + 8 <= q - p; i += 8) {
        __m512i hash = _mm512_setzero_si512();
        size_t o[8];
        size_t e[8];

        // Compute offsets and endpoints for each lane
        for (int lane = 0; lane < 8; ++lane) {
            if constexpr (level == 1)
                o[lane] = offsets1[p + i + lane];
            if constexpr (level == 2)
                o[lane] = offsets2[p + i + lane];
            if constexpr (level == 3)
                o[lane] = offsets3[p + i + lane];
            size_t next_endpoint = endpoints.select(endpoints.rank(o[lane]+1));
            size_t ee = o[lane] + k + span;
            e[lane] = (ee > next_endpoint) ? next_endpoint : ee;
        }

        // Compute initial hashes for each lane
        for (uint64_t j = 0; j < k; ++j) {
            uint64_t new_rank[8];
            for (int lane = 0; lane < 8; ++lane) {
                new_rank[lane] = seqan3::to_rank(text[o[lane] + j]);
            }
            __m512i ranks = _mm512_loadu_si512(new_rank);
            hash = _mm512_slli_epi64(hash, 2);
            hash = _mm512_or_si512(hash, ranks);
        }

        // Store initial hashes
        uint64_t hash_out[8];
        _mm512_storeu_si512(hash_out, hash);
        for (int lane = 0; lane < 8; ++lane) {
            buffer.push_back(hash_out[lane]);
        }

        // Continue for j = o[x]+k to e[x]
        __m512i mask_vec = _mm512_set1_epi64(mask);
        for (uint64_t j = k; j < k + span; ++j) {
            uint64_t new_rank[8];
            for (int lane = 0; lane < 8; ++lane) {
                if (o[lane] + j < e[lane])
                    new_rank[lane] = seqan3::to_rank(text[o[lane] + j]);
                else
                    new_rank[lane] = 0; // or some neutral value
            }
            __m512i ranks = _mm512_loadu_si512(new_rank);
            hash = _mm512_slli_epi64(hash, 2);
            hash = _mm512_or_si512(hash, ranks);
            hash = _mm512_and_si512(hash, mask_vec);

            _mm512_storeu_si512(hash_out, hash);
            for (int lane = 0; lane < 8; ++lane) {
                if (o[lane] + j < e[lane])
                    buffer.push_back(hash_out[lane]);
            }
        }
    }

    for (; i < q - p; ++i) {
        uint64_t hash = 0;
        size_t o;
        if constexpr (level == 1)
            o = offsets1[p + i];
        if constexpr (level == 2)
            o = offsets2[p + i];
        if constexpr (level == 3)
            o = offsets3[p + i];
        size_t next_endpoint = endpoints.select(endpoints.rank(o+1));
        size_t e = o + span + k;
        if (e > next_endpoint)
            e = next_endpoint;
        for (uint64_t j = o; j < o + k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
        }
        buffer.push_back(hash);
        for (size_t j = o + k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buffer.push_back(hash);
        }
    }

}


inline bool extend(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward) {
    if(forward) {
        if(last_found == array.size()-1)
            return false;
        if(array[last_found+1] == query) {
            last_found++;
            return true;
        }
    }
    else {
        if(last_found == 1)
            return false;
        if(array[last_found-1] == queryrc) {
            last_found--;
            return true;
        }
    }
    return false;
}


inline bool lookup_avx512(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward)
{
    const size_t n = array.size();

    __m512i qv = _mm512_set1_epi64(query);
    __m512i qrcv = _mm512_set1_epi64(queryrc);

    size_t i = 0;
    for (; i + 7 < n; i += 8) {
        __m512i v = _mm512_loadu_si512((const void*)&array[i]);
        __mmask8 mask1 = _mm512_cmpeq_epi64_mask(v, qv);
        __mmask8 mask2 = _mm512_cmpeq_epi64_mask(v, qrcv);
        if(mask1) {
            last_found = i + __builtin_clz(mask1);
            forward = true;
            return true;
        }
        if(mask2) {
            last_found = i + __builtin_clz(mask2);
            forward = false;
            return true;
        }
        // __mmask8 mask = mask1 | mask2;
        // if(mask) {
        //     int idx = __builtin_ctz(mask);
        //     forward = (mask1 & (1 << idx)) != 0;
        //     last_found = i + idx;
        //     return true;
        // }
    }
    for (; i < n; ++i) {
        if(array[i] == query) {
            last_found = i;
            forward = true;
            return true;
        }
        if(array[i] == queryrc) {
            last_found = i;
            forward = false;
            return true;
        }
    }

    return false;
}


inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward)
{
    for(size_t i=0; i < array.size(); i++) {
        if(array[i] == query) {
            last_found = i;
            forward = true;
            return true;
        }
        if(array[i] == queryrc) {
            last_found = i;
            forward = false;
            return true;
        }
    }
    return false;
}

inline bool lookup(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward, uint64_t &extensions)
{
    if(extend(array, query, queryrc, last_found, forward)) {
        extensions++;
        return true;
    }
    else
        // return lookup_serial(array, query, queryrc, last_found, forward);
        return lookup_avx512(array, query, queryrc, last_found, forward);
}


uint64_t RSIndexComp::streaming_query(const std::vector<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = srindex::views::three_minimisers_and_window_hash({.minimiser_size1 = m1, .minimiser_size2 = m2, .minimiser_size3 = m3, .window_size = k, .seed1=seed1, .seed2=seed2, .seed3=seed3});

    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(2u * k);
    uint64_t current_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser3=std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> buffer1;
    std::vector<uint64_t> buffer2;
    std::vector<uint64_t> buffer3;
    size_t last_found1 = 0;
    size_t last_found2 = 0;
    size_t last_found3 = 0;
    bool forward = true;

    for(auto && minimisers : query | view)
    {
        if(minimisers.minimiser1_value == current_minimiser1)
            occurences += lookup(buffer1, minimisers.window_value, minimisers.window_value_rev, last_found1, forward, extensions);
        else if(r1[minimisers.minimiser1_value]) {
            size_t minimizer_id = r1_rank(minimisers.minimiser1_value);
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);

            buffer1.clear();
            // fill_buffer<1>(buffer1, mask, p, q);
            fill_buffer_avx512<1>(buffer1, mask, p, q);
            last_found1 = 0;
            occurences += lookup(buffer1, minimisers.window_value, minimisers.window_value_rev, last_found1, forward, extensions);
            current_minimiser1 = minimisers.minimiser1_value;
        }
        else if(minimisers.minimiser2_value == current_minimiser2)
            occurences += lookup(buffer2, minimisers.window_value, minimisers.window_value_rev, last_found2, forward, extensions);
        else if(r2[minimisers.minimiser2_value]) {
            size_t minimizer_id = r2_rank(minimisers.minimiser2_value);
            size_t p = s2_select.select(minimizer_id);
            size_t q = s2_select.select(minimizer_id+1);

            buffer2.clear();
            // fill_buffer<2>(buffer2, mask, p, q);
            fill_buffer_avx512<2>(buffer2, mask, p, q);
            last_found2 = 0;
            occurences += lookup(buffer2, minimisers.window_value, minimisers.window_value_rev, last_found2, forward, extensions);
            current_minimiser2 = minimisers.minimiser2_value;
        }
        else if(minimisers.minimiser3_value == current_minimiser3)
            occurences += lookup(buffer3, minimisers.window_value, minimisers.window_value_rev, last_found3, forward, extensions);
        else if(r3[minimisers.minimiser3_value]) {
            size_t minimizer_id = r3_rank(minimisers.minimiser3_value);
            size_t p = s3_select.select(minimizer_id);
            size_t q = s3_select.select(minimizer_id+1);

            buffer3.clear();
            // fill_buffer<3>(buffer3, mask, p, q);
            fill_buffer_avx512<3>(buffer3, mask, p, q);
            last_found3 = 0;
            occurences += lookup(buffer3, minimisers.window_value, minimisers.window_value_rev, last_found3, forward, extensions);
            current_minimiser3 = minimisers.minimiser3_value;
        }
        else
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window_value, minimisers.window_value_rev));
    }
    
    return occurences;
}


int RSIndexComp::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m1, out);
    seqan3::contrib::sdsl::serialize(this->m2, out);
    seqan3::contrib::sdsl::serialize(this->m3, out);
    seqan3::contrib::sdsl::serialize(r1, out);
    seqan3::contrib::sdsl::serialize(r2, out);
    seqan3::contrib::sdsl::serialize(r3, out);
    seqan3::contrib::sdsl::serialize(s1, out);
    seqan3::contrib::sdsl::serialize(s2, out);
    seqan3::contrib::sdsl::serialize(s3, out);
    seqan3::contrib::sdsl::serialize(this->offsets1, out);
    seqan3::contrib::sdsl::serialize(this->offsets2, out);
    seqan3::contrib::sdsl::serialize(this->offsets3, out);
    seqan3::contrib::sdsl::serialize(this->sequences, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->text);
    archive(this->hashmap);

    out.close();
    return 0;
}

int RSIndexComp::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m1, in);
    seqan3::contrib::sdsl::load(this->m2, in);
    seqan3::contrib::sdsl::load(this->m3, in);
    seqan3::contrib::sdsl::load(r1, in);
    seqan3::contrib::sdsl::load(r2, in);
    seqan3::contrib::sdsl::load(r3, in);
    seqan3::contrib::sdsl::load(s1, in);
    seqan3::contrib::sdsl::load(s2, in);
    seqan3::contrib::sdsl::load(s3, in);
    seqan3::contrib::sdsl::load(this->offsets1, in);
    seqan3::contrib::sdsl::load(this->offsets2, in);
    seqan3::contrib::sdsl::load(this->offsets3, in);
    seqan3::contrib::sdsl::load(this->sequences, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);
    archive(this->hashmap);

    in.close();

    r1_rank = rank_support_v<1>(&r1);
    r2_rank = rank_support_sd<>(&r2);
    r3_rank = rank_support_sd<>(&r3);
    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);
    this->s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), s3.size(), 3);

    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), sequences.size());
    
    return 0;
}

