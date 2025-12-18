#include <filesystem>
#include <bitset>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
}


RSHash3::RSHash3() : endpoints(std::vector<uint64_t>{}, 1) {}

RSHash3::RSHash3(
    uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m3,
    uint8_t const m_thres1, uint8_t const m_thres2, uint16_t const m_thres3)
    : k(k), m1(m1), m2(m2), m3(m3),
      m_thres1(m_thres1), m_thres2(m_thres2), m_thres3(m_thres3), span1(k-m1+1), span2(k-m2+1), span3(k-m3+1),
      endpoints(std::vector<uint64_t>{}, 1)
{}


int RSHash3::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    auto view3 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m3, .window_size = k, .seed=seed3});
    auto view4 = srindex::views::xor_minimiser_and_window({.minimiser_size = m3, .window_size = k, .seed=seed3});

    std::cout << +m1 << " " << +m2 << " " << +m3 << "\n";

    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);
    const uint64_t M3 = 1ULL << (m3+m3);

    std::cout << "scan text...\n";
    size_t N = 0;
    uint64_t kmers = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        N += record.size();
        kmers += record.size() - k + 1;
        no_sequences++;
    }

    std::cout << "text length: " << N << "\n";
    std::cout << "text kmers: " << kmers <<  '\n';
    std::cout << "no sequences: " << no_sequences << "\n";

    std::cout << "get sequences...\n";
    sequences = bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), N+1);

    std::cout << "count minimizers...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers1;

    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            // minimizers1[minimiser.minimiser_value] += (minimiser.occurrences+span1-1)/span1;
            minimizers1[minimiser.minimiser_value]++;
            if(minimizers1[minimiser.minimiser_value] > m_thres1)
                minimizers1[minimiser.minimiser_value] = m_thres1;
            n++;
        }
    }
    
    std::cout << "extract unfrequent minimizers...\n";
    uint64_t n1 = 0;
    std::vector<uint64_t> unfreq_minimizers1;
    for(auto const& [minimizer, count] : minimizers1) {
        if(count < m_thres1) {
            unfreq_minimizers1.push_back(minimizer);
            n1 += count;
        }
    }

    size_t c1 = minimizers1.size();
    size_t c1tmp = unfreq_minimizers1.size();

    std::cout << "unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%)\n";

    std::cout << "build R_1...\n";
    sd_vector_builder builder(M1, unfreq_minimizers1.size());

    std::sort(unfreq_minimizers1.begin(), unfreq_minimizers1.end());
    for(uint64_t minimizer : unfreq_minimizers1)
        builder.set(minimizer);

    r1 = sd_vector<>(builder);
    r1_rank = rank_support_sd<>(&r1);

    std::cout << "filling bitvector S_1...\n";
    s1 = bit_vector(n1+1, 0);
    s1[0] = 1;
    j = 0;
    for(uint64_t minimizer : unfreq_minimizers1) {
        j += minimizers1[minimizer];
        s1[j] = 1;
    }
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), n1+1, 3);

    minimizers1.clear();
    unfreq_minimizers1.clear();

    std::cout << "filling offsets_1...\n";
    const size_t offset_width = std::bit_width(N);
    pthash::compact_vector::builder b1;
    b1.resize(n1, offset_width);

    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                size_t s = s1_select.select(i);
                b1.set(s + count1[i], length + minimiser.range_position);
                count1[i]++;
            }
        }
        length += sequence.size();
    }
    b1.build(offsets1);

    delete[] count1;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers1;
    std::vector<size_t> skmer_positions;
    length = 0;
    for(auto & sequence : input) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | view1) {
            freq = r1[minimiser.minimiser_value];
            break;
        }
        for(auto && minimiser : sequence | view1) {
            cur_freq = r1[minimiser.minimiser_value];
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers1.push_back(skmer);
                skmer_positions.push_back(length + start_position);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers1.push_back(skmer);
            skmer_positions.push_back(length + start_position);
        }
        length += sequence.size();
    }

    // uint64_t freq_kmers = 0;
    // for(auto & sequence : input) {
    //     for(auto && minimiser : sequence | srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1})) {
    //         freq_kmers += !r1[minimiser.minimiser_value];
    //     }
    // }
    // std::cout << "frequent k-mers: " << freq_kmers << " (" << (double) freq_kmers/kmers*100 << "%)\n";
    
    size_t rem_kmers1 = 0;
    for(auto & skmer : freq_skmers1)
        rem_kmers1 += skmer.size() - k + 1;
    std::cout << "remaining kmers: " << rem_kmers1 << " (" << (double) rem_kmers1/kmers*100 << "%)\n";

    std::cout << "build level 2...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers2;

    for(auto & skmer : freq_skmers1) {
        for(auto && minimiser : skmer | view2) {
            // minimizers2[minimiser.minimiser_value] += (minimiser.occurrences+span2-1)/span2;
            minimizers2[minimiser.minimiser_value]++;
            if(minimizers2[minimiser.minimiser_value] > m_thres2)
                minimizers2[minimiser.minimiser_value] = m_thres2;
        }
    }
    
    uint64_t n2 = 0;
    std::vector<uint64_t> unfreq_minimizers2;
    for(auto const& [minimizer, count] : minimizers2) {
        if(count < m_thres2) {
            unfreq_minimizers2.push_back(minimizer);
            n2 += count;
        }
    }

    size_t c2tmp = unfreq_minimizers2.size();

    std::cout << "build R_2...\n";
    std::sort(unfreq_minimizers2.begin(), unfreq_minimizers2.end());

    sd_vector_builder builder2(M2, unfreq_minimizers2.size());

    for(uint64_t minimizer : unfreq_minimizers2)
        builder2.set(minimizer);

    r2 = sd_vector<>(builder2);
    r2_rank = rank_support_sd<>(&r2);

    std::cout << "filling bitvector S_2...\n";
    s2 = bit_vector(n2+1, 0);
    s2[0] = 1;
    j = 0;
    for(size_t i=0; i < unfreq_minimizers2.size(); i++) {
        j += minimizers2[unfreq_minimizers2[i]];
        s2[j] = 1;
    }
    s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), n2+1, 3);

    minimizers2.clear();
    unfreq_minimizers2.clear();

    std::cout << "filling offsets_2...\n";
    pthash::compact_vector::builder b2;
    b2.resize(n2, offset_width);

    size_t c2 = r2_rank(M2);
    uint8_t* count2 = new uint8_t[c2];
    std::memset(count2, 0, c2*sizeof(uint8_t));

    size_t skmer_idx = 0;
    for(auto & skmer : freq_skmers1) {
        for(auto && minimiser : skmer | view2) {
            if(r2[minimiser.minimiser_value]) {
                size_t i = r2_rank(minimiser.minimiser_value);
                size_t s = s2_select.select(i);
                b2.set(s + count2[i], skmer_positions[skmer_idx] + minimiser.range_position);
                count2[i]++;
            }
        }
        skmer_idx++;
    }
    b2.build(offsets2);

    delete[] count2;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers2;
    std::vector<size_t> skmer_positions2;
    skmer_idx = 0; 
    for(auto & sequence : freq_skmers1) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | view2) {
            freq = r2[minimiser.minimiser_value];
            break;
        }
        for(auto && minimiser : sequence | view2) {
            cur_freq = r2[minimiser.minimiser_value];

            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers2.push_back(skmer);
                skmer_positions2.push_back(skmer_positions[skmer_idx] + start_position);
            }

            freq = cur_freq;
        }
        if(!cur_freq) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers2.push_back(skmer);
            skmer_positions2.push_back(skmer_positions[skmer_idx] + start_position);
        }

        skmer_idx++;
    }
    
    size_t rem_kmers2 = 0;
    for(auto & skmer : freq_skmers2)
        rem_kmers2 += skmer.size() - k + 1;
    // std::cout << "skmers: " << skmers << " minimizer: " << n << " unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%) \n";
    // std::cout << "remaining superkmers " << freq_skmers2.size() << " (" << (double) freq_skmers2.size()/skmers*100 << "%) "
    std::cout << "remaining kmers: " << rem_kmers2 << " (" << (double) rem_kmers2/kmers*100 << "%)\n";

    std::cout << "build level 3...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers3;

    for(auto & skmer : freq_skmers2) {
        for(auto && minimiser : skmer | view3) {
            // minimizers3[minimiser.minimiser_value] += (minimiser.occurrences+span3-1)/span3;
            minimizers3[minimiser.minimiser_value]++;
            if(minimizers3[minimiser.minimiser_value] > m_thres3)
                minimizers3[minimiser.minimiser_value] = m_thres3;
        }
    }
    
    uint64_t n3 = 0;
    std::vector<uint64_t> unfreq_minimizers3;
    for(auto const& [minimizer, count] : minimizers3) {
        if(count < m_thres3) {
            unfreq_minimizers3.push_back(minimizer);
            n3 += count;
        }
    }

    size_t c3tmp = unfreq_minimizers3.size();

    std::cout << "build R_3...\n";
    std::sort(unfreq_minimizers3.begin(), unfreq_minimizers3.end());

    sd_vector_builder builder3(M3, unfreq_minimizers3.size());

    for(uint64_t minimizer : unfreq_minimizers3)
        builder3.set(minimizer);

    r3 = sd_vector<>(builder3);
    r3_rank = rank_support_sd<>(&r3);

    std::cout << "filling bitvector S_3...\n";
    s3 = bit_vector(n3+1, 0);
    s3[0] = 1;
    j = 0;
    for(size_t i=0; i < unfreq_minimizers3.size(); i++) {
        j += minimizers3[unfreq_minimizers3[i]];
        s3[j] = 1;
    }
    s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), n3+1, 3);

    minimizers3.clear();
    unfreq_minimizers3.clear();

    std::cout << "filling offsets_3...\n";
    pthash::compact_vector::builder b3;
    b3.resize(n3, offset_width);

    size_t c3 = r3_rank(M3);
    uint8_t* count3 = new uint8_t[c3];
    std::memset(count3, 0, c3*sizeof(uint8_t));

    skmer_idx = 0;
    for(auto & skmer : freq_skmers2) {
        for(auto && minimiser : skmer | view3) {
            if(r3[minimiser.minimiser_value]) {
                size_t i = r3_rank(minimiser.minimiser_value);
                size_t s = s3_select.select(i);
                b3.set(s + count3[i], skmer_positions2[skmer_idx] + minimiser.range_position);
                count3[i]++;
            }
        }
        skmer_idx++;
    }
    b3.build(offsets3);

    delete[] count3;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers3;
    std::vector<size_t> skmer_positions3;
    skmer_idx = 0; 
    for(auto & sequence : freq_skmers2) {
        size_t start_position = 0;
        bool cur_freq, freq;

        for(auto && minimiser : sequence | view3) {
            freq = r3[minimiser.minimiser_value];
            break;
        }
        for(auto && minimiser : sequence | view3) {
            cur_freq = r3[minimiser.minimiser_value];
            if(freq && !cur_freq)
                start_position = minimiser.range_position;
            if(!freq && cur_freq) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position-1+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers3.push_back(skmer);
                skmer_positions3.push_back(skmer_positions2[skmer_idx] + start_position);
            }
            freq = cur_freq;
        }
        if(!cur_freq) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers3.push_back(skmer);
            skmer_positions3.push_back(skmer_positions2[skmer_idx] + start_position);
        }

        skmer_idx++;
    }

    std::cout << "build level 4, HT...\n";

    size_t rem_kmers3 = 0;
    for(auto & skmer : freq_skmers3)
        rem_kmers3 += skmer.size() - k + 1;
    // std::cout << "skmers: " << skmers << " minimizer: " << n << " unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%) \n";
    // std::cout << "remaining superkmers " << freq_skmers2.size() << " (" << (double) freq_skmers2.size()/skmers*100 << "%) "
    std::cout << "remaining kmers: " << rem_kmers3 << " (" << (double) rem_kmers3/kmers*100 << "%)\n";

    std::cout << "filling HT...\n";
    // todo: simple kmer view
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
    std::cout << "Hashtable mem: " << (double) 65*hashmap.bucket_count()/kmers << "\n";
    std::cout << "R_1: " << (double) 8*size_in_bytes(r1)/kmers << "\n";
    std::cout << "R_2: " << (double) 8*size_in_bytes(r2)/kmers << "\n";
    std::cout << "R_3: " << (double) 8*size_in_bytes(r3)/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1)/kmers << "\n";
    std::cout << "S_2: " << (double) (n2+1)/kmers << "\n";
    std::cout << "S_3: " << (double) (n3+1)/kmers << "\n";

    std::cout << "total: " << (double) (n1*offset_width+n2*offset_width+n3*offset_width+2*N+8*size_in_bytes(r1)+8*size_in_bytes(r2)+8*size_in_bytes(r3)+n1+1+n2+1+n3+1+endpoints.bitCount()+64*hashmap.size())/kmers << "\n";
    std::cout << "total mem: " << (double) (n1*offset_width+n2*offset_width+n3*offset_width+2*N+8*size_in_bytes(r1)+8*size_in_bytes(r2)+8*size_in_bytes(r3)+n1+1+n2+1+n3+1+endpoints.bitCount()+65*hashmap.bucket_count())/kmers << "\n";

    return 0;
}



static inline constexpr uint64_t crc(uint64_t x, uint64_t k) {
    // assert(k <= 32);
    uint64_t c = ~x;

    /* swap byte order */
    uint64_t res = __builtin_bswap64(c);

    /* Swap nuc order in bytes */
    const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;              // ...0000.1111.0000.1111
    const uint64_t c2 = 0x3333333333333333;              // ...0011.0011.0011.0011
    res = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4);  // swap 2-nuc order in bytes
    res = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2);  // swap nuc order in 2-nuc

    /* Realign to the right */
    res >>= 64 - 2 * k;

    return res;
}


std::vector<uint64_t> RSHash3::rand_text_kmers(const uint64_t n) {
    std::uniform_int_distribution<uint32_t> distr;
    std::mt19937 m_rand(1);
    std::vector<std::uint64_t> kmers;
    kmers.reserve(n);

    const uint64_t no_unitigs = number_unitigs();
    for (uint64_t i = 0; i < n; ++i) {
        const uint64_t unitig_id = distr(m_rand) % no_unitigs;
        const uint64_t offset = distr(m_rand) % unitig_size(unitig_id);
        const uint64_t kmer = access(unitig_id, offset);

        if ((i & 1) == 0)
            kmers.push_back(crc(kmer, k));
        else
            kmers.push_back(kmer);
    }

    return kmers;
}



uint64_t RSHash3::access(const uint64_t unitig_id, const size_t offset)
{
    size_t offset_text = endpoints.select(unitig_id) + offset;

    uint64_t kmer = 0;
    for (size_t i=offset_text; i < offset_text+k; i++) {
        uint64_t const new_rank = seqan3::to_rank(text[i]);
        kmer <<= 2;
        kmer |= new_rank;
    }

    return kmer;
}


template<int level>
inline bool RSHash3::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc,
    double &to, double &th, double &te)
{
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    for(size_t i = 0; i < q-p; i++) {

        t0 = std::chrono::high_resolution_clock::now();
        size_t o, span;
        if constexpr (level == 1) {
            o = offsets1.access(p+i);
            span = span1;
        }
        if constexpr (level == 2) {
            o = offsets2.access(p+i);
            span = span2;
        }
        if constexpr (level == 3) {
            o = offsets3.access(p+i);
            span = span3;
        }
        t1 = std::chrono::high_resolution_clock::now();
        to += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();

        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        t2 = std::chrono::high_resolution_clock::now();
        th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
        if(hash == kmer || hash == kmer_rc)
            return true;

        uint64_t e = std::min<uint64_t>(endpoints.select(endpoints.rank(o+1)), o+k+span);
        // uint64_t e = o+k+span;

        t3 = std::chrono::high_resolution_clock::now();
        te += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();

        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            if(hash == kmer || hash == kmer_rc) {
                t4 = std::chrono::high_resolution_clock::now();
                th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
                return true;
            }
        }
        t4 = std::chrono::high_resolution_clock::now();
        th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
    }

    return false;
}


template<int level>
inline bool RSHash3::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc)
{
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    for(size_t i = 0; i < q-p; i++) {

        t0 = std::chrono::high_resolution_clock::now();
        size_t o, span;
        if constexpr (level == 1) {
            o = offsets1.access(p+i);
            span = span1;
        }
        if constexpr (level == 2) {
            o = offsets2.access(p+i);
            span = span2;
        }
        if constexpr (level == 3) {
            o = offsets3.access(p+i);
            span = span3;
        }

        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        if(hash == kmer || hash == kmer_rc)
            return true;

        uint64_t e = std::min<uint64_t>(endpoints.select(endpoints.rank(o+1)), o+k+span);
        // uint64_t e = o+k+span;

        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            if(hash == kmer || hash == kmer_rc) {
                return true;
            }
        }
    }

    return false;
}

uint64_t RSHash3::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
uint64_t occurences = 0;
if (verbose) {
    const uint64_t mask = compute_mask(2u * k);
    srindex::minimizers::Three_minimisers_hash2 minimisers = srindex::minimizers::Three_minimisers_hash2(k, m1, m2, m3, seed1, seed2, seed3);
    std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4, t5, t6 = std::chrono::high_resolution_clock::now();
    double t0_ = 0.0;
    double t1_ = 0.0;
    double t2_ = 0.0;
    double t3_ = 0.0;
    double t4_ = 0.0;
    double t5_ = 0.0;
    double to = 0.0;
    double th = 0.0;
    double te = 0.0;
    size_t skmers1_ = 0;
    size_t skmers2_ = 0;
    size_t skmers3_ = 0;
    uint64_t lookups1 = 0;
    uint64_t lookups2 = 0;
    uint64_t lookups3 = 0;
    uint64_t ht_lookups = 0;

    for(uint64_t kmer : kmers)
    {
        t5 = std::chrono::high_resolution_clock::now();
        minimisers.compute(kmer);
        t4_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t5)).count();

        t6 = std::chrono::high_resolution_clock::now();
        if(uint64_t minimizer_id = r1_rank(minimisers.minimiser1); r1_rank(minimisers.minimiser1+1) > minimizer_id) {
        // if(r1[minimisers.minimiser1]) {
            t0 = std::chrono::high_resolution_clock::now();
            // uint64_t minimizer_id = r1_rank(minimisers.minimiser1);
            t1 = std::chrono::high_resolution_clock::now();
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);
            t2 = std::chrono::high_resolution_clock::now();

            occurences += check<1>(p, q, mask, minimisers.window, minimisers.window_rev, to, th, te);
            // occurences += check<1>(p, q, mask, minimisers.window, minimisers.window_rev);
            t3 = std::chrono::high_resolution_clock::now();

            t0_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();
            t1_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            t2_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - t6)).count();
            skmers1_ += q - p;
            ++lookups1;
        }
        else if(r2[minimisers.minimiser2]) {
            t0 = std::chrono::high_resolution_clock::now();
            uint64_t minimizer_id = r2_rank(minimisers.minimiser2);
            t1 = std::chrono::high_resolution_clock::now();
            size_t p = s2_select.select(minimizer_id);
            size_t q = s2_select.select(minimizer_id+1);
            t2 = std::chrono::high_resolution_clock::now();

            occurences += check<2>(p, q, mask, minimisers.window, minimisers.window_rev, to, th, te);
            // occurences += check<2>(p, q, mask, minimisers.window, minimisers.window_rev);
            t3 = std::chrono::high_resolution_clock::now();

            t0_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();
            t1_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            t2_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - t6)).count();
            skmers2_ += q - p;
            ++lookups2;
        }
        else if(r3[minimisers.minimiser3]) {
            t0 = std::chrono::high_resolution_clock::now();
            uint64_t minimizer_id = r3_rank(minimisers.minimiser3);
            t1 = std::chrono::high_resolution_clock::now();
            size_t p = s3_select.select(minimizer_id);
            size_t q = s3_select.select(minimizer_id+1);
            t2 = std::chrono::high_resolution_clock::now();

            occurences += check<3>(p, q, mask, minimisers.window, minimisers.window_rev, to, th, te);
            // occurences += check<3>(p, q, mask, minimisers.window, minimisers.window_rev);
            t3 = std::chrono::high_resolution_clock::now();

            t0_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();
            t1_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            t2_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - t6)).count();
            skmers3_ += q - p;
            ++lookups3;
        }
        else {
            t4 = std::chrono::high_resolution_clock::now();
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window, minimisers.window_rev));
            t3_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t4)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t6)).count();
            ++ht_lookups;
        }
    }
    std::cout << "minimisers: " << t4_/kmers.size() << " ns\n";
    std::cout << "R lookup: " << t5_/kmers.size() << " ns\n";
    std::cout << "r_rank: " << t0_/kmers.size() << " ns\n";
    std::cout << "s_select: " << t1_/kmers.size() << " ns\n";
    std::cout << "check: " << t2_/kmers.size() << " ns\n";
    std::cout << "offsets: " << to/kmers.size() << " ns\n";
    std::cout << "endpoints: " << te/kmers.size() << " ns\n";
    std::cout << "text: " << th/kmers.size() << " ns\n";
    std::cout << "ht: " << t3_/kmers.size() << " ns\n";
    std::cout << "lookups lvl1: " << (double) lookups1/(lookups1+lookups2+lookups3+ht_lookups)*100 << "%\n";
    std::cout << "lookups lvl2: " << (double) lookups2/(lookups1+lookups2+lookups3+ht_lookups)*100 << "%\n";
    std::cout << "lookups lvl3: " << (double) lookups3/(lookups1+lookups2+lookups3+ht_lookups)*100 << "%\n";
    std::cout << "lookups ht: " << (double) ht_lookups/(lookups1+lookups2+lookups3+ht_lookups)*100 << "%\n";
    std::cout << "avg skmers1: " << (double) skmers1_/lookups1 << "\n";
    std::cout << "avg skmers2: " << (double) skmers2_/lookups2 << "\n";
    std::cout << "avg skmers3: " << (double) skmers3_/lookups3 << "\n";
}
else {
    const uint64_t mask = compute_mask(2u * k);
    srindex::minimizers::Three_minimisers_hash2 minimisers = srindex::minimizers::Three_minimisers_hash2(k, m1, m2, m3, seed1, seed2, seed3);

    for(uint64_t kmer : kmers)
    {
        minimisers.compute(kmer);

        if(r1[minimisers.minimiser1]) {
            uint64_t minimizer_id = r1_rank(minimisers.minimiser1);
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);

            occurences += check<1>(p, q, mask, minimisers.window, minimisers.window_rev);
        }
        else if(r2[minimisers.minimiser2]) {
            uint64_t minimizer_id = r2_rank(minimisers.minimiser2);
            size_t p = s2_select.select(minimizer_id);
            size_t q = s2_select.select(minimizer_id+1);

            occurences += check<2>(p, q, mask, minimisers.window, minimisers.window_rev);
        }
        else if(r3[minimisers.minimiser3]) {
            uint64_t minimizer_id = r3_rank(minimisers.minimiser3);
            size_t p = s3_select.select(minimizer_id);
            size_t q = s3_select.select(minimizer_id+1);

            occurences += check<3>(p, q, mask, minimisers.window, minimisers.window_rev);
        }
        else {
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window, minimisers.window_rev));
        }

    }
}

    return occurences;
}





template<int level>
inline void RSHash3::fill_buffer(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers, size_t p, size_t q)
{
    for(uint64_t i = 0; i < q-p; i++) {
        uint64_t hash = 0;
        size_t o, span;
        if constexpr (level == 1) {
            o = offsets1.access(p+i);
            span = span1;
        }
        if constexpr (level == 2) {
            o = offsets2.access(p+i);
            span = span2;
        }
        if constexpr (level == 3) {
            o = offsets3.access(p+i);
            span = span3;
        }
        uint64_t next_endpoint;
        uint64_t prev_endpoint = endpoints.select(endpoints.rank(o+1)-1, &next_endpoint);
        size_t e = std::min(o+k+span, next_endpoint);
        
        if(e > next_endpoint)
            e = next_endpoint;
        for (uint64_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        buffer.push_back(hash);
        for (uint64_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            buffer.push_back(hash);
        }
        skmers.push_back({o, prev_endpoint, next_endpoint});
    }
}



inline bool RSHash3::lookup_serial(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers,
    const uint64_t query, const uint64_t queryrc,
    size_t &text_pos, bool &forward, size_t &start_pos, size_t &end_pos)
{
    size_t s = 0, e = 0;
    for (size_t i = 0; i < skmers.size(); i++) {
        e += 1;
        for(size_t j = s; j < e; j++) {
            if(buffer[j] == query) {
                forward = true;
                text_pos = skmers[i].position + (j - s) + k - 1;
                end_pos = skmers[i].unitig_end;
                return true;
            }
            if(buffer[j] == queryrc) {
                forward = false;
                text_pos = skmers[i].position + (j - s);
                start_pos = skmers[i].unitig_begin;
                return true;
            }
        }
        s = e;
    }
    return false;
}


inline bool RSHash3::extend_in_text(size_t &text_pos, size_t start, size_t end,
    bool forward, const uint64_t query, const uint64_t query_rc)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query >> (2*(k-1))));
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query_rc & 0b11));
            return found;
        }
    }
    return false;
}


uint64_t RSHash3::streaming_query(const std::vector<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = srindex::views::xor_three_minimiser_and_window({.minimiser1_size = m1, .minimiser2_size = m2, .minimiser3_size = m3, .window_size = k, .seed1=seed1, .seed2=seed2, .seed3=seed3});

    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(2u * k);
    uint64_t current_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser2=std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser3=std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> buffer1;
    std::vector<uint64_t> buffer2;
    std::vector<uint64_t> buffer3;
    std::vector<SkmerInfo> skmers1;
    std::vector<SkmerInfo> skmers2;
    std::vector<SkmerInfo> skmers3;
    size_t unitig_begin, unitig_end;
    size_t text_pos;
    bool forward;
    bool found = false;

    for(auto && minimisers : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, minimisers.window_value, minimisers.window_value_rev)) {
            occurences++;
            extensions++;
        }
        else if(minimisers.minimiser1_value == current_minimiser1) {
            found = lookup_serial(buffer1, skmers1, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(r1[minimisers.minimiser1_value]) {
            uint64_t minimizer_id = r1_rank(minimisers.minimiser1_value);
            size_t p = s1_select.select(minimizer_id);
            size_t q = s1_select.select(minimizer_id+1);

            buffer1.clear();
            skmers1.clear();
            fill_buffer<1>(buffer1, skmers1, p, q);
            found = lookup_serial(buffer1, skmers1, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser1 = minimisers.minimiser1_value;
        }
        else if(minimisers.minimiser2_value == current_minimiser2) {
            found = lookup_serial(buffer2, skmers2, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(r2[minimisers.minimiser2_value]) {
            uint64_t minimizer_id = r2_rank(minimisers.minimiser2_value);
            size_t p = s2_select.select(minimizer_id);
            size_t q = s2_select.select(minimizer_id+1);

            buffer2.clear();
            skmers2.clear();
            fill_buffer<2>(buffer2, skmers2, p, q);
            found = lookup_serial(buffer2, skmers2, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser2 = minimisers.minimiser2_value;
        }
        else if(minimisers.minimiser3_value == current_minimiser3) {
            found = lookup_serial(buffer3, skmers3, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(r3[minimisers.minimiser3_value]) {
            uint64_t minimizer_id = r3_rank(minimisers.minimiser3_value);
            size_t p = s3_select.select(minimizer_id);
            size_t q = s3_select.select(minimizer_id+1);

            buffer3.clear();
            skmers3.clear();
            fill_buffer<3>(buffer3, skmers3, p, q);
            found = lookup_serial(buffer3, skmers3, minimisers.window_value, minimisers.window_value_rev, text_pos, forward, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser3 = minimisers.minimiser3_value;
        }
        else {
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window_value, minimisers.window_value_rev));
            found = false;
        }
    }
    
    return occurences;
}


int RSHash3::save(const std::filesystem::path &filepath) {
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
    seqan3::contrib::sdsl::serialize(this->sequences, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->offsets1);
    archive(this->offsets2);
    archive(this->offsets3);
    archive(this->text);
    archive(this->hashmap);

    out.close();
    return 0;
}

int RSHash3::load(const std::filesystem::path &filepath) {
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
    seqan3::contrib::sdsl::load(this->sequences, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->offsets1);
    archive(this->offsets2);
    archive(this->offsets3);
    archive(this->text);
    archive(this->hashmap);

    this->span1 = k-m1+1;
    this->span2 = k-m2+1;
    this->span3 = k-m3+1;

    std::cout << "loaded index...\n";

    in.close();

    r1_rank = rank_support_sd<>(&r1);
    r2_rank = rank_support_sd<>(&r2);
    r3_rank = rank_support_sd<>(&r3);
    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);
    this->s3_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s3.data()), s3.size(), 3);

    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), sequences.size());

    std::cout << "built rank and select ds...\n";
    
    return 0;
}

