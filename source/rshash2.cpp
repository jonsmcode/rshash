#include <filesystem>
#include <bitset>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <cereal/archives/binary.hpp>

#include "rshash.hpp"
#include "io.hpp"
#include "minimiser_views.hpp"


RSHash2::RSHash2() :
    endpoints(std::vector<uint64_t>{}, 1),
    r1(std::vector<uint64_t>{}, 1),
    r2(std::vector<uint64_t>{}, 1)
{}

RSHash2::RSHash2(
    uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m_thres1, uint8_t const m_thres2)
    : k(k), m1(m1), m_thres1(m_thres1), m2(m2),
      m_thres2(m_thres2), span1(k-m1+1), span2(k-m2+1),
      endpoints(std::vector<uint64_t>{}, 1),
      r1(std::vector<uint64_t>{}, 1),
      r2(std::vector<uint64_t>{}, 1)
{}


int RSHash2::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto minimiserview1 = srindex::views::xor_minimiser_and_positions2({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto minimiserview2 = srindex::views::xor_minimiser_and_positions2({.minimiser_size = m2, .window_size = k, .seed=seed2});
    auto skmerview1 = srindex::views::xor_minimiser_and_skmer_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto skmerview2 = srindex::views::xor_minimiser_and_skmer_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    auto kmerview = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});

    const uint64_t M1 = 1ULL << (m1+m1);
    const uint64_t M2 = 1ULL << (m2+m2);

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
    bit_vector sequences = bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), N+1);

    std::cout << "count minimizers1...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers1;

    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | minimiserview1) {
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
    const size_t offset_width = std::bit_width(N);
    pthash::compact_vector::builder b1;
    b1.resize(n1, offset_width);

    uint8_t* count1 = new uint8_t[c1tmp];
    std::memset(count1, 0, c1tmp*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | minimiserview1) {
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

        for(auto && minimiser : sequence | skmerview1) {
            freq = r1.rank(minimiser.minimiser_value+1)-r1.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | skmerview1) {
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


    std::cout << "count minimizers2...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers2;

    for(auto & skmer : freq_skmers1) {
        for(auto && minimiser : skmer | minimiserview2) {
            minimizers2[minimiser.minimiser_value]++;
            if(minimizers2[minimiser.minimiser_value] > m_thres2)
                minimizers2[minimiser.minimiser_value] = m_thres2;
        }
    }
    
    std::cout << "extract unfrequent minimizers...\n";
    uint64_t n2 = 0;
    std::vector<uint64_t> unfreq_minimizers2;
    for(auto const& [minimizer, count] : minimizers2) {
        if(count < m_thres2) {
            unfreq_minimizers2.push_back(minimizer);
            n2 += count;
        }
    }

    size_t c2 = minimizers2.size();
    size_t c2tmp = unfreq_minimizers2.size();

    std::cout << "unfrequent minimizers: " << unfreq_minimizers2.size() << " (" << (double) unfreq_minimizers2.size()/minimizers2.size()*100 << "%)\n";

    std::cout << "build R_2...\n";
    std::sort(unfreq_minimizers2.begin(), unfreq_minimizers2.end());
    r2 = sux::bits::EliasFano(unfreq_minimizers2, M2);

    std::cout << "filling bitvector S_2...\n";
    s2 = bit_vector(n2+1, 0);
    s2[0] = 1;
    j = 0;
    for(uint64_t minimizer : unfreq_minimizers2) {
        j += minimizers2[minimizer];
        s2[j] = 1;
    }
    s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), n2+1, 3);

    minimizers2.clear();
    unfreq_minimizers2.clear();

    std::cout << "filling offsets_2...\n";
    pthash::compact_vector::builder b2;
    b2.resize(n2, offset_width);

    uint8_t* count2 = new uint8_t[c2tmp];
    std::memset(count2, 0, c2tmp*sizeof(uint8_t));

    size_t skmer_idx = 0;
    for(auto & skmer : freq_skmers1) {
        for (auto && minimiser : skmer | minimiserview2) {
            if(uint64_t i = r2.rank(minimiser.minimiser_value); r2.rank(minimiser.minimiser_value+1)-i) {
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

        for(auto && minimiser : sequence | skmerview2) {
            freq = r2.rank(minimiser.minimiser_value+1)-r2.rank(minimiser.minimiser_value);
            break;
        }
        for(auto && minimiser : sequence | skmerview2) {
            cur_freq = r2.rank(minimiser.minimiser_value+1)-r2.rank(minimiser.minimiser_value);
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
    // std::cout << "remaining superkmers " << freq_skmers1.size() << " (" << (double) freq_skmers1.size()/skmers*100 << "%) ";
    std::cout << "remaining kmers: " << rem_kmers2 << " (" << (double) rem_kmers2/kmers*100 << "%)\n";


    std::cout << "build level 3, HT...\n";

    // todo: simple kmer view
    // todo: hashmap.reserve(rem_kmers2);
    for(auto & sequence : freq_skmers2) {
        for(auto && minimiser : sequence | kmerview) {
            hashmap.insert(std::min<uint64_t>(minimiser.window_value, minimiser.window_value_rev));
        }
    }

    std::cout << "copy text...\n";
    text = pack_dna4_to_uint64(input);

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "textkmers: " << kmers <<  '\n';
    
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << c1tmp << "\n";
    std::cout << "minimiser going to level 2: " << c1-c1tmp << "  " << (double) (c1-c1tmp)/c1tmp*100 << "%\n";
    std::cout << "no minimiser1: " << n1 << "\n";
    std::cout << "no distinct minimiser1: " << c1 << "\n";
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    std::cout << "no minimiser2: " << n2 << "\n";
    std::cout << "no distinct minimiser2: " << c2 << "\n";
    std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/kmers*100 << "%\n";

    std::cout << "density r1: " << (double) c1/M1*100 << "%\n";
    std::cout << "density r2: " << (double) c2/M2*100 << "%\n";
    std::cout << "density s1: " << (double) s1_select.bitCount()/(n1+1)*100 <<  "%\n";
    std::cout << "density s2: " << (double) s2_select.bitCount()/(n2+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    std::cout << "endpoints: " << (double) endpoints.bitCount()/kmers << "\n";
    std::cout << "offsets1: " << (double) n1*offset_width/kmers << "\n";
    std::cout << "offsets2: " << (double) n2*offset_width/kmers << "\n";
    std::cout << "Hashtable: " << (double) 65*hashmap.bucket_count()/kmers << "\n";
    std::cout << "R_1: " << (double) r1.bitCount()/kmers << "\n";
    std::cout << "R_2: " << (double) r2.bitCount()/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1)/kmers << "\n";
    std::cout << "S_2: " << (double) (n2+1)/kmers << "\n";
    
    std::cout << "total: " << (double) (n1*offset_width+n2*offset_width+2*N+r1.bitCount()+r2.bitCount()+n1+1+n2+1+endpoints.bitCount()+65*hashmap.bucket_count())/kmers << "\n";

    return 0;
}


std::vector<uint64_t> RSHash2::rand_text_kmers(const uint64_t n) {
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



// uint64_t RSHash2::access(const uint64_t unitig_id, const size_t offset)
// {
//     uint64_t offset_text = endpoints.select(unitig_id) + offset;
//     const uint64_t mask = compute_mask(2u * k);
//     return get_word64(offset_text) & mask;
// }
uint64_t RSHash2::access(const uint64_t unitig_id, const size_t offset)
{
    const uint64_t mask = compute_mask(2u * k);
    return get_word64(offset) & mask;
}

template<int level>
inline bool RSHash2::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc,
    double &to, double &th, double &te)
{
    return false;
}

template<int level>
inline bool RSHash2::check(uint64_t* offsets, std::array<uint64_t, 2>* sequence_ends, const size_t p, const size_t no_kmers, const uint64_t mask, const uint64_t shift,
    const uint64_t kmer, const uint64_t kmer_rc)
{
    uint64_t span;
    if constexpr (level == 1)
        span = span1;
    if constexpr (level == 2)
        span = span2;

    for(size_t i = 0; i < no_kmers; i++) {
        if constexpr (level == 1)
            offsets[i] = offsets1.access(p+i);
        if constexpr (level == 2)
            offsets[i] = offsets2.access(p+i);
    }
    
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

uint64_t RSHash2::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
    uint64_t occurences = 0;

    const uint64_t mask = compute_mask(2u * k);
    const uint64_t shift = 2*(k-1);
    srindex::minimizers::Three_minimisers_hash2 minimisers = srindex::minimizers::Three_minimisers_hash2(k, m1, m2, 0, seed1, seed2, seed3);
    uint64_t* offsets = new uint64_t[m_thres2];
    std::array<uint64_t, 2>* sequence_ends = new std::array<uint64_t, 2>[m_thres2];

    for(uint64_t kmer : kmers)
    {
        const uint64_t kmer_rc = crc(kmer, k);
        minimisers.compute(kmer, kmer_rc);

        if(uint64_t minimiser_rank = r1.rank(minimisers.minimiser1); r1.rank(minimisers.minimiser1 + 1) - minimiser_rank) {
            size_t p = s1_select.select(minimiser_rank);
            size_t no_kmers = s1_select.select(minimiser_rank+1) - p;

            occurences += check<1>(offsets, sequence_ends, p, no_kmers, mask, shift, kmer, kmer_rc);
        }
        else if(uint64_t minimiser_rank = r2.rank(minimisers.minimiser2); r2.rank(minimisers.minimiser2 + 1) - minimiser_rank) {
            size_t p = s2_select.select(minimiser_rank);
            size_t no_kmers = s2_select.select(minimiser_rank+1) - p;

            occurences += check<2>(offsets, sequence_ends, p, no_kmers, mask, shift, kmer, kmer_rc);
        }
        else
            occurences += hashmap.contains(std::min<uint64_t>(kmer, kmer_rc));
    }

    delete[] offsets;
    delete[] sequence_ends;

    return occurences;
}


const inline uint64_t RSHash2::get_word64(uint64_t pos) {
    uint64_t block = pos >> 5;
    uint64_t shift = (pos & 31) << 1;
    uint64_t lo = text[block];
    uint64_t hi = text[block + 1];

    uint64_t shift_mask = -(shift != 0);
    return (lo >> shift) | ((hi << (64 - shift)) & shift_mask);
}

const inline uint64_t RSHash2::get_base(uint64_t pos) {
    return (text[pos >> 5] >> ((pos & 31) << 1)) & 3ULL;
}



inline bool RSHash2::extend_in_text(size_t &text_pos, size_t start, size_t end,
    bool forward, const uint64_t query, const uint64_t query_rc)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query >> (2*(k-1))));
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = get_base(text_pos);
            bool const found = (new_rank == (query_rc & 0b11));
            return found;
        }
    }
    return false;
}


template<int level>
inline void RSHash2::refill_buffer(uint64_t *buffer, SkmerInfo *skmers, size_t p, size_t N, const uint64_t mask, const uint64_t shift)
{
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();

    for(size_t i = 0; i < N; i++) {
        uint64_t o, span;
        if constexpr (level == 1) {
            o = offsets1.access(p+i);
            span = span1;
        }
        if constexpr (level == 2) {
            o = offsets2.access(p+i);
            span = span2;
        }

        uint64_t next_endpoint;
        const uint64_t r = endpoints.rank(o+1);
        const uint64_t prev_endpoint = endpoints.select(r-1, &next_endpoint);

        const uint64_t delta = std::min<uint64_t>(o+1, span);
        const uint64_t s1 = o + 1 - delta;
        const uint64_t s2 = std::max<uint64_t>(s1, prev_endpoint);
        const uint64_t e = std::min<uint64_t>(o+k, next_endpoint);

        skmers[i] = {s1, prev_endpoint, next_endpoint};
        
        const uint64_t pad_front = span - delta;
        std::fill(buffer, buffer + pad_front, INF);
        buffer += pad_front;
        if(prev_endpoint > s1) {
            const uint64_t pad_front2 = prev_endpoint - s1;
            std::fill(buffer, buffer + pad_front2, INF);
            buffer += pad_front2;
        }
        
        uint64_t kmer = get_word64(s2) & mask;
        uint64_t bits = get_word64(s2 + k);
        *buffer++ = kmer;
        for(uint64_t j=s2+k; j < e; j++) {
            uint64_t const next_base = bits & 3ULL;
            bits >>= 2;
            kmer = (kmer >> 2) | (next_base << shift);
            *buffer++ = kmer;
        }

        const uint64_t pad_back = o + k - e;
        std::fill(buffer, buffer + pad_back, INF);
        buffer += pad_back;
    }
}

inline bool RSHash2::check_minimiser_pos(uint64_t* buffer, const SkmerInfo &skmer,
    const uint64_t query, const uint64_t queryrc,
    const size_t s, const size_t e, const size_t minimiser_pos,
    bool &forward, size_t &text_pos, size_t &start_pos, size_t &end_pos)
{
    if(buffer[s+minimiser_pos] == query) {
        forward = true;
        text_pos = skmer.position + minimiser_pos + k - 1;
        end_pos = skmer.unitig_end;
        return true;
    }
    if(buffer[s+minimiser_pos] == queryrc) {
        forward = false;
        text_pos = skmer.position + minimiser_pos;
        start_pos = skmer.unitig_begin;
        return true;
    }
    if(buffer[e-1-minimiser_pos] == queryrc) {
        forward = false;
        text_pos = skmer.position + e-1-s-minimiser_pos;
        start_pos = skmer.unitig_begin;
        return true;
    }
    if(buffer[e-1-minimiser_pos] == query) {
        forward = true;
        text_pos = skmer.position + e-1-s-minimiser_pos + k - 1;
        end_pos = skmer.unitig_end;
        return true;
    }

    return false;
}


template<int level>
inline bool RSHash2::lookup_buffer(uint64_t* buffer, SkmerInfo* skmers, const size_t no_skmers,
    const uint64_t query, const uint64_t queryrc,
    size_t &text_pos, const size_t left_minimiser_pos, const size_t right_minimiser_pos,
    bool &forward, size_t &start_pos, size_t &end_pos)
{
    size_t span, m;
    if constexpr (level == 1) {
        span = span1;
        m = m1;
    }
    if constexpr (level == 2) {
        span = span2;
        m = m2;
    }

    size_t s = 0, e = 0;
    if(left_minimiser_pos != k-m-right_minimiser_pos) {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;
            if(check_minimiser_pos(buffer, skmers[i], query, queryrc, s, e, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            if(check_minimiser_pos(buffer, skmers[i], query, queryrc, s, e, right_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    else {
        for(size_t i = 0; i < no_skmers; i++) {
            e += span;   
            if(check_minimiser_pos(buffer, skmers[i], query, queryrc, s, e, left_minimiser_pos, forward, text_pos, start_pos, end_pos))
                return true;
            s = e;
        }
    }
    
    return false;
}


uint64_t RSHash2::streaming_query(const std::vector<seqan3::dna4> &query, uint64_t &extensions)
{
    auto view = srindex::views::xor_two_minimiser_and_window2({.minimiser1_size = m1, .minimiser2_size = m2, .window_size = k, .seed1=seed1, .seed2=seed2});

    const uint64_t mask = compute_mask(2u * k);
    const uint64_t shift = 2*(k-1);
    constexpr uint64_t INF = std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser1=INF, current_minimiser2=INF;
    uint64_t current_neg_minimiser1=INF, current_neg_minimiser2=INF;
    uint64_t* buffer1 = new uint64_t[(m_thres1-1) * span1];
    SkmerInfo* skmers1 = new SkmerInfo[m_thres1-1];
    uint64_t* buffer2 = new uint64_t[(m_thres2-1) * span2];
    SkmerInfo* skmers2 = new SkmerInfo[m_thres2-1];
    uint64_t minimiser_rank1, minimiser_rank2;
    size_t no_skmers1, no_skmers2;
    size_t unitig_begin, unitig_end;
    size_t text_pos;
    bool forward;
    bool found = false;
    uint64_t occurences = 0;

    for(auto && window : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, window.window_value, window.window_value_rev)) {
            occurences++;
            extensions++;
        }
        else if(window.minimiser1_value == current_minimiser1) {
            found = lookup_buffer<1>(buffer1, skmers1, no_skmers1, window.window_value, window.window_value_rev, text_pos, window.minimiser1_left_position, window.minimiser1_right_position, forward, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(window.minimiser1_value != current_neg_minimiser1 && (minimiser_rank1 = r1.rank(window.minimiser1_value), r1.rank(window.minimiser1_value + 1) - minimiser_rank1)) {
            const size_t p = s1_select.select(minimiser_rank1);
            no_skmers1 = s1_select.select(minimiser_rank1+1) - p;

            refill_buffer<1>(buffer1, skmers1, p, no_skmers1, mask, shift);
            found = lookup_buffer<1>(buffer1, skmers1, no_skmers1, window.window_value, window.window_value_rev, text_pos, window.minimiser1_left_position, window.minimiser1_right_position, forward, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser1 = window.minimiser1_value;
        }
        else if(window.minimiser2_value == current_minimiser2) {
            found = lookup_buffer<2>(buffer2, skmers2, no_skmers2, window.window_value, window.window_value_rev, text_pos, window.minimiser2_left_position, window.minimiser2_right_position, forward, unitig_begin, unitig_end);
            occurences += found;
        }
        else if(window.minimiser2_value != current_neg_minimiser2 && (minimiser_rank2 = r2.rank(window.minimiser2_value), r2.rank(window.minimiser2_value + 1) - minimiser_rank2)) {
            const size_t p = s2_select.select(minimiser_rank2);
            no_skmers2 = s2_select.select(minimiser_rank2+1) - p;

            refill_buffer<2>(buffer2, skmers2, p, no_skmers2, mask, shift);
            found = lookup_buffer<2>(buffer2, skmers2, no_skmers2, window.window_value, window.window_value_rev, text_pos, window.minimiser2_left_position, window.minimiser2_right_position, forward, unitig_begin, unitig_end);
            occurences += found;
            current_minimiser2 = window.minimiser2_value;
            current_neg_minimiser1 = window.minimiser1_value;
        }
        else {
            occurences += hashmap.contains(std::min<uint64_t>(window.window_value, window.window_value_rev));
            found = false;
            current_neg_minimiser1 = window.minimiser1_value;
            current_neg_minimiser2 = window.minimiser2_value;
        }

    }

    delete[] buffer1;
    delete[] buffer2;
    delete[] skmers1;
    delete[] skmers2;
    
    return occurences;
}




int RSHash2::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m1, out);
    seqan3::contrib::sdsl::serialize(this->m2, out);
    seqan3::contrib::sdsl::serialize(this->m_thres1, out);
    seqan3::contrib::sdsl::serialize(this->m_thres2, out);
    seqan3::contrib::sdsl::serialize(s1, out);
    seqan3::contrib::sdsl::serialize(s2, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->endpoints);
    archive(this->r1);
    archive(this->r2);
    archive(this->offsets1);
    archive(this->offsets2);
    archive(this->text);
    archive(this->hashmap);

    out.close();
    return 0;
}

int RSHash2::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m1, in);
    seqan3::contrib::sdsl::load(this->m2, in);
    seqan3::contrib::sdsl::load(this->m_thres1, in);
    seqan3::contrib::sdsl::load(this->m_thres2, in);
    seqan3::contrib::sdsl::load(s1, in);
    seqan3::contrib::sdsl::load(s2, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->endpoints);
    archive(this->r1);
    archive(this->r2);
    archive(this->offsets1);
    archive(this->offsets2);
    archive(this->text);
    archive(this->hashmap);

    this->span1 = k - m1 + 1;
    this->span2 = k - m2 + 1;

    std::cout << "loaded index...\n";

    in.close();

    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);

    std::cout << "built rank and select ds...\n";
    
    return 0;
}
