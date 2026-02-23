#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


seqan3::dna4_vector kmer_to_string(uint64_t kmer, size_t const kmer_size)
{
    seqan3::dna4_vector result(kmer_size);
    for (size_t i = 0; i < kmer_size; ++i)
    {
        result[kmer_size - 1 - i].assign_rank(kmer & 0b11);
        kmer >>= 2;
    }
    return result;
}


RSHash1::RSHash1() :
    endpoints(std::vector<uint64_t>{}, 1),
    r1(std::vector<uint64_t>{}, 1),
    m_hasher(seed1)
{}

RSHash1::RSHash1(
    uint8_t const k, uint8_t const m1, uint8_t const m_thres1)
    : k(k), m1(m1),
      m_thres1(m_thres1), span(k-m1+1),
      endpoints(std::vector<uint64_t>{}, 1),
      r1(std::vector<uint64_t>{}, 1),
      m_hasher(seed1)
{}


int RSHash1::build(std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions2({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view3 = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});

    std::cout << +m1 << "\n";

    const uint64_t M1 = 1ULL << (m1+m1);

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
    bit_vector sequences = bit_vector(N+2, 0);
    sequences[0] = 1;
    sequences[32] = 1;
    size_t j = 32;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), N+2);
    sequences = bit_vector();

    std::cout << "count minimizers1...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers1;

    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
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
    std::sort(unfreq_minimizers1.begin(), unfreq_minimizers1.end());

    std::cout << "maximum minimizer1 value: " << unfreq_minimizers1.back() << " max value " << M1 <<  " fraction " << (double) unfreq_minimizers1.back()/M1*100 << "%\n";

    r1 = sux::bits::EliasFano(unfreq_minimizers1, M1);

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
    const size_t offset_width = std::bit_width(N+span);
    pthash::compact_vector::builder b1;
    b1.resize(n1, offset_width);

    uint8_t* count1 = new uint8_t[c1tmp];
    std::memset(count1, 0, c1tmp*sizeof(uint8_t));

    size_t length = 32;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
            if(uint64_t i = r1.rank(minimiser.minimiser_value); r1.rank(minimiser.minimiser_value+1)-i) {
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

        for(auto && minimiser : sequence | view2) {
            freq = r1.rank(minimiser.minimiser_value+1)-r1.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | view2) {
            cur_freq = r1.rank(minimiser.minimiser_value+1)-r1.rank(minimiser.minimiser_value);
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
    
    size_t rem_kmers1 = 0;
    for(auto & skmer : freq_skmers1)
        rem_kmers1 += skmer.size() - k + 1;
    std::cout << "remaining kmers: " << rem_kmers1 << " (" << (double) rem_kmers1/kmers*100 << "%)\n";

    std::cout << "build level 2...\n";

    // size_t i = 0;
    // for(const auto& sequence : freq_skmers1) {
    //     std::cout << ">freq superkmer " << i++ << std::endl;
    //     for(const auto& nucleotide : sequence)
    //         std::cout << nucleotide.to_char();
    //     std::cout << std::endl;
    // }

    for(auto & sequence : freq_skmers1) {
        for(auto && window : sequence | srindex::views::kmerview({.window_size = k})) {
            hashmap.insert(std::min<uint64_t>(window.kmer_value, window.kmer_value_rev));
        }
    }

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "textkmers: " << kmers <<  '\n';
    
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << c1tmp << "\n";
    std::cout << "minimiser going to level 2: " << c1-c1tmp << "  " << (double) (c1-c1tmp)/c1*100 << "%\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1tmp << "\n";
    std::cout << "avg superkmers1: " << (double) n1/c1tmp <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/kmers*100 << "%\n";

    std::cout << "density r1: " << (double) c1/M1*100 << "%\n";
    std::cout << "density s1: " << (double) s1_select.bitCount()/(n1+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/kmers << "\n";
    std::cout << "offsets1: " << (double) n1*offset_width/kmers << "\n";
    std::cout << "R_1: " << (double) r1.bitCount()/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1 + s1_select.bitCount())/kmers << "\n";
    std::cout << "Hashtable: " << (double) 65*hashmap.bucket_count()/kmers << "\n";

    std::cout << "total: " << (double) (n1*offset_width+2*N+r1.bitCount()+n1+1+s1_select.bitCount()+endpoints.bitCount()+65*hashmap.bucket_count())/kmers << "\n";

    return 0;
}


std::vector<uint64_t> RSHash1::rand_text_kmers(const uint64_t n) {
    std::uniform_int_distribution<uint32_t> distr;
    std::mt19937 m_rand(1);
    std::vector<std::uint64_t> kmers;
    kmers.reserve(n);
    const uint64_t l = (text.size()-1)*32;

    const uint64_t no_unitigs = number_unitigs();
    for (uint64_t i = 0; i < n;) {
        const uint64_t offset = distr(m_rand) % l;

        const uint64_t r = endpoints.rank(offset+1);
        const uint64_t next_endpoint = endpoints.select(r);

        if(offset + 64 >= next_endpoint)
            continue;
        // const uint64_t unitig_id = distr(m_rand) % no_unitigs;
        // const uint64_t offset = distr(m_rand) % unitig_size(unitig_id);
        const uint64_t kmer = access(0, offset);

        if ((i & 1) == 0)
            kmers.push_back(crc(kmer, k));
        else
            kmers.push_back(kmer);

        i++;
    }

    return kmers;
}



uint64_t RSHash1::access(const uint64_t unitig_id, const size_t offset)
{
    const uint64_t mask = compute_mask(2u * k);
    return get_word64(offset) & mask;
}


inline bool RSHash1::check(const size_t p, const size_t no_kmers, const uint64_t mask, const uint64_t shift,
    const uint64_t kmer, const uint64_t kmer_rc,
    double &to, double &th, double &te)
{
    // std::chrono::high_resolution_clock::time_point t0, t1, t2, t3, t4;
    // for(size_t i = 0; i < q-p; i++) {

    //     t0 = std::chrono::high_resolution_clock::now();
    //     size_t o = offsets1.access(p+i);
    //     t1 = std::chrono::high_resolution_clock::now();
    //     to += (std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0)).count();

    //     uint64_t hash = 0;
    //     for (size_t j=o; j < o+k; j++) {
    //         uint64_t const new_rank = seqan3::to_rank(text[j]);
    //         hash = (hash >> 2) | (new_rank << 2*(k-1));
    //     }
    //     t2 = std::chrono::high_resolution_clock::now();
    //     th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
    //     if(hash == kmer || hash == kmer_rc)
    //         return true;

    //     uint64_t e = std::min<uint64_t>(endpoints.select(endpoints.rank(o+1)), o+k+span);
    //     // uint64_t e = o+k+span;

    //     t3 = std::chrono::high_resolution_clock::now();
    //     te += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();

    //     for(size_t j=o+k; j < e; j++) {
    //         uint64_t const new_rank = seqan3::to_rank(text[j]);
    //         hash = (hash >> 2) | (new_rank << 2*(k-1));
    //         if(hash == kmer || hash == kmer_rc) {
    //             t4 = std::chrono::high_resolution_clock::now();
    //             th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
    //             return true;
    //         }
    //     }
    //     t4 = std::chrono::high_resolution_clock::now();
    //     th += (std::chrono::duration_cast<std::chrono::nanoseconds>(t4 - t3)).count();
    // }

    return false;
}


inline bool RSHash1::check(uint64_t* offsets, std::array<uint64_t, 2>* sequence_ends, const size_t p, const size_t no_kmers, const uint64_t mask, const uint64_t shift,
    const uint64_t kmer, const uint64_t kmer_rc)
{
    for(size_t i = 0; i < no_kmers; i++)
        offsets[i] = offsets1.access(p+i);
    for(size_t i = 0; i < no_kmers; i++) {
        const uint64_t o = offsets[i];
        uint64_t next_endpoint;
        const uint64_t r = endpoints.rank(o+1);
        const uint64_t prev_endpoint = endpoints.select(r-1, &next_endpoint);

        const uint64_t delta = std::min<uint64_t>(o+1, span);
        const uint64_t s1 = o + 1 - delta;
        const uint64_t s2 = std::max<uint64_t>(s1, prev_endpoint);
        const uint64_t e = std::min<uint64_t>(o+k, next_endpoint);

        sequence_ends[i] = {s2, e};
    }
    for(size_t i = 0; i < no_kmers; i++) {
        auto [s2, e] = sequence_ends[i];
        uint64_t hash = get_word64(s2) & mask;
        if(hash == kmer || hash == kmer_rc)
            return true;

        uint64_t bits = get_word64(s2 + k);
        for(uint64_t j=s2+k; j < e; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            hash = (hash >> 2) | (next_base << shift);
            if(hash == kmer || hash == kmer_rc)
                return true;
        }
    }

    return false;
}

uint64_t RSHash1::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(2u * k);
    const uint64_t shift = 2*(k-1);
    srindex::minimizers::Minimisers_hash2 minimisers = srindex::minimizers::Minimisers_hash2(k, m1, seed1);
    uint64_t* offsets = new uint64_t[m_thres1];
    std::array<uint64_t, 2>* sequence_ends = new std::array<uint64_t, 2>[m_thres1];

    if (verbose) {
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
        uint64_t lookups1 = 0;
        uint64_t ht_lookups = 0;

        for(uint64_t kmer : kmers)
        {
        t5 = std::chrono::high_resolution_clock::now();
        const uint64_t kmer_rc = crc(kmer, k);
        uint64_t minimiser_value = minimisers.compute(kmer, kmer_rc);
        t4_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - t5)).count();

        t6 = std::chrono::high_resolution_clock::now();
        if(uint64_t minimiser_rank = r1.rank(minimiser_value); r1.rank(minimiser_value + 1) - minimiser_rank) {
            t1 = std::chrono::high_resolution_clock::now();
            size_t p = s1_select.select(minimiser_rank);
            size_t no_kmers = s1_select.select(minimiser_rank+1) - p;
            t2 = std::chrono::high_resolution_clock::now();

            occurences += check(p, no_kmers, mask, shift, kmer, kmer_rc, to, th, te);
            t3 = std::chrono::high_resolution_clock::now();

            t1_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1)).count();
            t2_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2)).count();
            t5_ += (std::chrono::duration_cast<std::chrono::nanoseconds>(t0 - t6)).count();
            skmers1_ += no_kmers;
            ++lookups1;
        }
        else {
            t4 = std::chrono::high_resolution_clock::now();
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
        std::cout << "lookups lvl1: " << (double) lookups1/(lookups1+ht_lookups)*100 << "%\n";
        std::cout << "lookups ht: " << (double) ht_lookups/(lookups1+ht_lookups)*100 << "%\n";
        std::cout << "avg skmers1: " << (double) skmers1_/lookups1 << "\n";
    }
    else {
        for(uint64_t kmer : kmers)
        {
            const uint64_t kmer_rc = crc(kmer, k);
            uint64_t minimiser_value = minimisers.compute(kmer, kmer_rc);

            if(uint64_t minimiser_rank = r1.rank(minimiser_value); r1.rank(minimiser_value + 1) - minimiser_rank) {
                size_t p = s1_select.select(minimiser_rank);
                size_t no_kmers = s1_select.select(minimiser_rank+1) - p;

                occurences += check(offsets, sequence_ends, p, no_kmers, mask, shift, kmer, kmer_rc);
            }
            else
                occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
        }
    }

    delete[] offsets;
    delete[] sequence_ends;

    return occurences;
}




inline bool RSHash1::extend_in_text(size_t &text_pos, size_t start, size_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, const uint64_t shift)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = get_base(text_pos);
            return new_rank == (query >> shift);
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = get_base(text_pos);
            return new_rank == (query_rc & 0b11);
        }
    }
    return false;
}


const inline uint64_t RSHash1::get_word64(uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t RSHash1::get_base(uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
}


inline void RSHash1::fill_buffer(uint64_t *offsets, uint64_t *buffer, size_t p, size_t N, const uint64_t mask, const uint64_t shift)
{
    for(uint64_t i = 0; i < N; i++)
        offsets[i] = offsets1.access(p+i)+1-span;

    for(uint64_t i = 0; i < N; i++) {
        const uint64_t o = offsets[i];
        
        uint64_t kmer = get_word64(o) & mask;
        uint64_t bits = get_word64(o + k);
        *buffer++ = kmer;
        for(size_t j=0; j < span-1; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            kmer = (kmer >> 2) | (next_base << shift);
            *buffer++ = kmer;
        }
    }
}

inline bool RSHash1::check_overlap(uint64_t skmer_pos, uint64_t text_pos, uint64_t &start_pos, uint64_t &end_pos)
{
    const uint64_t r = endpoints.rank(skmer_pos+span);
    start_pos = endpoints.select(r-1, &end_pos);

    return text_pos >= start_pos && text_pos+k-1 < end_pos;
}


inline bool RSHash1::check_minimiser_pos(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, const size_t s,
    const size_t minimiser_pos, bool &forward, size_t &text_pos, size_t &start_pos, size_t &end_pos)
{
    if(buffer[s+minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + span-1-minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }

    return false;
}

inline bool RSHash1::check_minimiser_pos2(uint64_t *buffer, const uint64_t offset,
    const uint64_t kmer, const uint64_t kmerrc, const size_t s,
    const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, size_t &text_pos, size_t &start_pos, size_t &end_pos)
{
    if(buffer[s+left_minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + left_minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-left_minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + span-1-left_minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+right_minimiser_pos] == kmer) {
        forward = true;
        text_pos = offset + right_minimiser_pos + k - 1;
        if(check_overlap(offset, text_pos-k+1, start_pos, end_pos))
            return true;
    }
    if(buffer[s+span-1-right_minimiser_pos] == kmerrc) {
        forward = false;
        text_pos = offset + span-1-right_minimiser_pos;
        if(check_overlap(offset, text_pos, start_pos, end_pos))
            return true;
    }

    return false;
}


inline bool RSHash1::lookup_buffer(uint64_t *buffer, uint64_t *offsets, const size_t no_minimiser,
    const uint64_t query, const uint64_t queryrc,
    size_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, size_t &start_pos, size_t &end_pos)
{
    size_t s = 0;
    if(left_minimiser_pos != k-m1-right_minimiser_pos) {
        for(size_t i = 0; i < no_minimiser; i++) {
            if(check_minimiser_pos2(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, right_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s += span;
        }
    }
    else {
        for(size_t i = 0; i < no_minimiser; i++) {
            if(check_minimiser_pos(buffer, offsets[i], query, queryrc, s, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s += span;
        }
    }
    
    return false;
}


inline uint64_t RSHash1::find_minimiser(const uint64_t kmer, const uint64_t kmer_rc, size_t &left_minimiser_position, size_t &right_minimiser_position, const uint64_t mmermask)
{
    uint64_t mmer = kmer >> 2*(k - m1);
    uint64_t mmer_rc = kmer_rc & mmermask;
    uint64_t minimiser_value = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
    left_minimiser_position = k-m1;
    right_minimiser_position = 0;

    for (size_t i = 1; i <= k-m1; ++i) {
        mmer = (kmer >> 2*(k-m1-i)) & mmermask;
        mmer_rc = (kmer_rc >> 2*i) & mmermask;
        const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);
        if(mmerhash < minimiser_value) {
            minimiser_value = mmerhash;
            left_minimiser_position = k-m1-i;
            right_minimiser_position = i;
        }
        else if(mmerhash == minimiser_value)
            left_minimiser_position = k-m1-i;
    }

    return minimiser_value;
}

inline void RSHash1::update_minimiser(const uint64_t kmer, const uint64_t kmer_rc, uint64_t &minimiser, size_t &left_minimiser_position, size_t &right_minimiser_position, const uint64_t mmermask)
{
    if(left_minimiser_position-- == 0) {
        minimiser = find_minimiser(kmer, kmer_rc, left_minimiser_position, right_minimiser_position, mmermask);
        return;
    }
    const uint64_t mmer = kmer >> 2*(k - m1);
    const uint64_t mmer_rc = kmer_rc & mmermask;
    const uint64_t mmerhash = std::min<uint64_t>(m_hasher.hash(mmer) & mmermask, m_hasher.hash(mmer_rc) & mmermask);

    if(mmerhash < minimiser) {
        minimiser = mmerhash;
        left_minimiser_position = k - m1;
        right_minimiser_position = 0;
        return;
    }
    if(mmerhash == minimiser) {
        right_minimiser_position = 0;
        return;
    }

    right_minimiser_position++;
}


uint64_t RSHash1::streaming_query(const seqan3::bitpacked_sequence<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = srindex::views::kmerview({.window_size = k});

    uint64_t occurences = 0;
    uint64_t current_pos_minimiser=std::numeric_limits<uint64_t>::max();
    uint64_t current_neg_minimiser=std::numeric_limits<uint64_t>::max();
    const uint64_t kmermask = compute_mask(2u * k);
    const uint64_t mmermask = compute_mask(2u * m1);
    const uint64_t shift = 2*(k-1);
    uint64_t* offsets = new uint64_t[m_thres1-1];
    uint64_t* kmer_buffer = new uint64_t[(m_thres1-1) * span];
    size_t no_minimiser, text_pos, sequence_begin, sequence_end;
    bool forward;
    bool found = false;
    bool rolling = false;
    uint64_t minimiser, minimiser_rank;
    size_t left_minimiser_position, right_minimiser_position;

    for(auto && kmer : query | view)
    {
        if(found && extend_in_text(text_pos, sequence_begin, sequence_end, forward, kmer.kmer_value, kmer.kmer_value_rev, shift)) {
            occurences++;
            extensions++;
            rolling = false;
        }
        else {
            if(rolling)
                update_minimiser(kmer.kmer_value, kmer.kmer_value_rev, minimiser, left_minimiser_position, right_minimiser_position, mmermask);
            else {
                minimiser = find_minimiser(kmer.kmer_value, kmer.kmer_value_rev, left_minimiser_position, right_minimiser_position, mmermask);
                rolling = true;
            }

            if(minimiser == current_pos_minimiser) {
                found = lookup_buffer(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, text_pos, left_minimiser_position, right_minimiser_position, forward, sequence_begin, sequence_end);
                occurences += found;
            }
            else if(minimiser != current_neg_minimiser && (minimiser_rank = r1.rank(minimiser), r1.rank(minimiser + 1) - minimiser_rank)) {
                const size_t minimiser_position = s1_select.select(minimiser_rank);
                no_minimiser = s1_select.select(minimiser_rank+1) - minimiser_position;

                fill_buffer(offsets, kmer_buffer, minimiser_position, no_minimiser, kmermask, shift);
                found = lookup_buffer(kmer_buffer, offsets, no_minimiser, kmer.kmer_value, kmer.kmer_value_rev, text_pos, left_minimiser_position, right_minimiser_position, forward, sequence_begin, sequence_end);
                occurences += found;
                current_pos_minimiser = minimiser;
            }
            else {
                occurences += hashmap.contains(std::min<uint64_t>(kmer.kmer_value, kmer.kmer_value_rev));
                found = false;
                current_neg_minimiser = minimiser;
            }
        }

    }

    delete[] kmer_buffer;
    delete[] offsets;
    
    return occurences;
}


void RSHash1::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    cereal::BinaryOutputArchive archive(out);

    archive(k, m1, m_thres1, s1, endpoints, r1, offsets1, text, hashmap);

    out.close();
}

void RSHash1::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    cereal::BinaryInputArchive archive(in);

    archive(k, m1, m_thres1, s1, endpoints, r1, offsets1, text, hashmap);

    in.close();

    span = k - m1 + 1;
    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);

    std::cout << "loaded index...\n";
}

