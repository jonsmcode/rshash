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


RSHash2::RSHash2() : endpoints(std::vector<uint64_t>{}, 1) {}

RSHash2::RSHash2(
    uint8_t const k, uint8_t const m1, uint8_t const m2, uint8_t const m_thres1, uint8_t const m_thres2)
    : k(k), m1(m1), m_thres1(m_thres1), m2(m2),
      m_thres2(m_thres2), span1(k-m1+2), span2(k-m2+2),
      endpoints(std::vector<uint64_t>{}, 1)
{}


int RSHash2::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view2 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m2, .window_size = k, .seed=seed2});
    auto view3 = srindex::views::xor_minimiser_and_window({.minimiser_size = m1, .window_size = k, .seed=seed1});
    auto view4 = srindex::views::xor_minimiser_and_window({.minimiser_size = m2, .window_size = k, .seed=seed2});

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
    sequences = bit_vector(N+1, 0);
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
        for(auto && minimiser : sequence | view1) {
            minimizers1[minimiser.minimiser_value] += minimiser.occurrences/span1+1;
            if(minimizers1[minimiser.minimiser_value] > m_thres1)
                minimizers1[minimiser.minimiser_value] = m_thres1;
            n += minimiser.occurrences/span1+1;
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

    uint8_t* count1 = new uint8_t[c1tmp];
    std::memset(count1, 0, c1tmp*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                size_t s = s1_select.select(i);
                size_t o = minimiser.occurrences;
                size_t j = 0;
                while(o > span1) {
                    b1.set(s + count1[i], length + minimiser.range_position + j*span1);
                    count1[i]++;
                    o -= span1;
                    j++;
                }
                b1.set(s + count1[i], length + minimiser.range_position + j*span1);
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

    uint64_t freq_kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view3) {
            freq_kmers += !r1[minimiser.minimiser_value];
        }
    }
    std::cout << "frequent k-mers: " << freq_kmers << " (" << (double) freq_kmers/kmers*100 << "%)\n";
    
    size_t rem_kmers1 = 0;
    for(auto & skmer : freq_skmers1)
        rem_kmers1 += skmer.size() - k + 1;
    // std::cout << "skmers: " << skmers << " minimizer: " << n << " unfrequent minimizers: " << unfreq_minimizers1.size() << " (" << (double) unfreq_minimizers1.size()/minimizers1.size()*100 << "%) \n";
    // std::cout << "remaining superkmers " << freq_skmers1.size() << " (" << (double) freq_skmers1.size()/skmers*100 << "%) ";
    std::cout << "remaining kmers: " << rem_kmers1 << " (" << (double) rem_kmers1/kmers*100 << "%)\n";


    std::cout << "count minimizers2...\n";
    std::unordered_map<uint64_t, uint8_t> minimizers2;

    for(auto & skmer : freq_skmers1) {
        for(auto && minimiser : skmer | view2) {
            minimizers2[minimiser.minimiser_value] += minimiser.occurrences/span2+1;
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
    sd_vector_builder builder2(M2, unfreq_minimizers2.size());

    std::sort(unfreq_minimizers2.begin(), unfreq_minimizers2.end());
    for(uint64_t minimizer : unfreq_minimizers2)
        builder2.set(minimizer);

    r2 = sd_vector<>(builder2);
    r2_rank = rank_support_sd<>(&r2);

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
        for (auto && minimiser : skmer | view2) {
            if(r2[minimiser.minimiser_value]) {
                size_t i = r2_rank(minimiser.minimiser_value);
                size_t s = s2_select.select(i);
                size_t o = minimiser.occurrences;
                size_t j = 0;
                while(o > span2) {
                    b2.set(s + count2[i], skmer_positions[skmer_idx] + minimiser.range_position + j*span2);
                    count2[i]++;
                    o -= span2;
                    j++;
                }
                b2.set(s + count2[i], skmer_positions[skmer_idx] + minimiser.range_position + j*span2);
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

    freq_kmers = 0;
    for(auto & sequence : freq_skmers1) {
        for(auto && minimiser : sequence | view4) {
            freq_kmers += !r2[minimiser.minimiser_value];
        }
    }
    std::cout << "frequent k-mers: " << freq_kmers << " (" << (double) freq_kmers/kmers*100 << "%)\n";
    
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
        for(auto && minimiser : sequence | view3) {
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
    std::cout << "Hashtable: " << (double) 64*hashmap.size()/kmers << "\n";
    // std::cout << "Hashtable: " << (double) 64*hashmap.bucket_count()/kmers << "\n";
    std::cout << "Hashtable mem: " << (double) 65*hashmap.bucket_count()/kmers << "\n";
    std::cout << "R_1: " << (double) 8*size_in_bytes(r1)/kmers << "\n";
    std::cout << "R_2: " << (double) 8*size_in_bytes(r2)/kmers << "\n";
    std::cout << "S_1: " << (double) (n1+1)/kmers << "\n";
    std::cout << "S_2: " << (double) (n2+1)/kmers << "\n";
    
    // todo: bitvektor rank select overhead
    std::cout << "total: " << (double) (n1*offset_width+n2*offset_width+2*N+8*size_in_bytes(r1)+8*size_in_bytes(r2)+n1+1+n2+1+endpoints.bitCount()+64*hashmap.size())/kmers << "\n";

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


std::vector<uint64_t> RSHash2::rand_text_kmers(const uint64_t n) {
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



uint64_t RSHash2::access(const uint64_t unitig_id, const size_t offset)
{
    size_t offset_text = endpoints.select(unitig_id) + offset;

    uint64_t kmer = 0;
    // for (size_t i=offset_text; i < offset_text+k; i++) {
    //     uint64_t const new_rank = seqan3::to_rank(text[i]);
    //     kmer <<= 2;
    //     kmer |= new_rank;
    // }

    return kmer;
}

template<int level>
inline bool RSHash2::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc,
    double &to, double &th, double &te)
{
    return false;
}

template<int level>
inline bool RSHash2::check(const size_t p, const size_t q, const uint64_t mask,
    const uint64_t kmer, const uint64_t kmer_rc)
{
    return false;
}

uint64_t RSHash2::lookup(const std::vector<uint64_t> &kmers, bool verbose)
{
    uint64_t occurences = 0;
    return occurences;
}



// inline bool extend(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool forward) {
//     if(forward) {
//         // if(last_found == array.size()-1)
//         //     return false;
//         if(array[last_found+1] == query) {
//             last_found++;
//             return true;
//         }
//     }
//     else {
//         // if(last_found == 1)
//         //     return false;
//         if(array[last_found-1] == queryrc) {
//             last_found--;
//             return true;
//         }
//     }
//     return false;
// }


inline bool RSHash2::extend_in_text(size_t &text_pos, size_t start, size_t end,
    bool forward, const uint64_t query, const uint64_t query_rc, uint64_t &fwd_extensions, uint64_t &rev_extensions)
{
    if(forward) {
        if(++text_pos < end) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query >> (2*(k-1))));
            fwd_extensions++;
            return found;
        }
    }
    else {
        if(--text_pos >= start) {
            uint64_t const new_rank = seqan3::to_rank(text[text_pos]);
            bool const found = (new_rank == (query_rc & 0b11));
            rev_extensions++;
            return found;
        }
    }
    return false;
}


inline bool extend_in_buffer(std::vector<uint64_t> &buffer, const uint64_t query, const uint64_t queryrc,
    size_t &skmer_pos, bool forward, uint64_t &fwd_extensions, uint64_t &rev_extensions)
{
    if(forward) {
        if(buffer[skmer_pos+1] == query) {
            skmer_pos++;
            fwd_extensions++;
            return true;
        }
    }
    else {
        if(buffer[skmer_pos-1] == queryrc) {
            skmer_pos--;
            rev_extensions++;
            return true;
        }
    }
    return false;
}


template<int level>
inline void RSHash2::fill_buffer(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers, size_t p, size_t q)
{
    for(size_t i = 0; i < q-p; i++) {
        size_t o;
        size_t span;
        if constexpr (level == 1) {
            o = offsets1.access(p+i);
            span = span1;
        }
        if constexpr (level == 2) {
            o = offsets2.access(p+i);
            span = span2;
        }
        uint64_t next_endpoint;
        uint64_t prev_endpoint = endpoints.select(endpoints.rank(o+1)-1, &next_endpoint);
        size_t e = std::min(o+k+span, next_endpoint);

        // todo: get 256 chars around o in one cache line
        // todo: let buffer be bits text[max(o-126,s),...,o,...,min(o+126,e)]
        // check kmer at position p in skmer with reinterpret_cast<uint64_t*>(buffer+2*p)[0] >> (64 - 2*k);
        uint64_t hash = 0;
        for (size_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
        }
        buffer.push_back(hash);
        for (size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash = (hash >> 2) | (new_rank << 2*(k-1));
            buffer.push_back(hash);
        }
        skmers.push_back({o, e - o - k + 1, prev_endpoint, next_endpoint});
    }
}


inline bool RSHash2::lookup_serial(std::vector<uint64_t> &buffer, std::vector<SkmerInfo> &skmers,
    const uint64_t query, const uint64_t queryrc,
    size_t &text_pos, bool &forward, size_t &start_pos, size_t &end_pos)
{
    size_t s = 0, e = 0;
    for (size_t i = 0; i < skmers.size(); i++) {
        e += skmers[i].length;
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


// inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward)
// {
//     const size_t n = array.size();
//     for (int i = 0; i < n; ++i) {
//         if(array[i] == query) {
//             last_found = i;
//             forward = true;
//             return true;
//         }
//         if(array[i] == queryrc) {
//             last_found = i;
//             forward = false;
//             return true;
//         }
//     }
//     return false;
// }


// inline bool streaming_lookup(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found, bool &forward, uint64_t &extensions)
// {
//     if(extend(array, query, queryrc, last_found, forward)) {
//         extensions++;
//         return true;
//     }
//     else {
//         return lookup_serial(array, query, queryrc, last_found, forward);
//     }
        
// }

uint64_t RSHash2::streaming_query(const std::vector<seqan3::dna4> &query,
    uint64_t &buffer_fwd_extensions, uint64_t &buffer_rev_extensions, uint64_t &text_fwd_extensions, uint64_t &text_rev_extensions)
{
    auto view = srindex::views::xor_two_minimiser_and_window({.minimiser1_size = m1, .minimiser2_size = m2, .window_size = k, .seed1=seed1, .seed2=seed2});

    uint64_t occurences = 0;
    uint64_t current_minimiser1=std::numeric_limits<uint64_t>::max();
    uint64_t current_minimiser2=std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> buffer1;
    std::vector<uint64_t> buffer2;
    std::vector<SkmerInfo> skmers1;
    std::vector<SkmerInfo> skmers2;
    size_t unitig_begin, unitig_end;
    size_t text_pos;
    bool forward;
    bool found = false;

    for(auto && minimisers : query | view)
    {
        if(found && extend_in_text(text_pos, unitig_begin, unitig_end, forward, minimisers.window_value, minimisers.window_value_rev, text_fwd_extensions, text_rev_extensions))
            occurences++;
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
        else {
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window_value, minimisers.window_value_rev));
            found = false;
        }

    }
    
    return occurences;
}




int RSHash2::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m1, out);
    seqan3::contrib::sdsl::serialize(this->m2, out);
    seqan3::contrib::sdsl::serialize(r1, out);
    seqan3::contrib::sdsl::serialize(r2, out);
    seqan3::contrib::sdsl::serialize(s1, out);
    seqan3::contrib::sdsl::serialize(s2, out);
    seqan3::contrib::sdsl::serialize(this->sequences, out);

    cereal::BinaryOutputArchive archive(out);
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
    seqan3::contrib::sdsl::load(r1, in);
    seqan3::contrib::sdsl::load(r2, in);
    seqan3::contrib::sdsl::load(s1, in);
    seqan3::contrib::sdsl::load(s2, in);
    seqan3::contrib::sdsl::load(this->sequences, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->offsets1);
    archive(this->offsets2);
    archive(this->text);
    archive(this->hashmap);

    this->span1 = k - m1 + 2;
    this->span2 = k - m2 + 2;

    std::cout << "loaded index...\n";

    in.close();

    r1_rank = rank_support_sd<>(&r1);
    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    r2_rank = rank_support_sd<>(&r2);
    this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);

    endpoints = sux::bits::EliasFano(reinterpret_cast<uint64_t*>(sequences.data()), sequences.size());

    std::cout << "built rank and select ds...\n";
    
    return 0;
}
