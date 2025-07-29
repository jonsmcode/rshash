#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <cereal/archives/binary.hpp>

#include "dict.hpp"
#include "../source/minimiser_rev_xor_hash_views3.hpp"


const uint8_t m_thres1 = 10;
const uint8_t m_thres2 = 10;

const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 1;

const size_t span = 100;


static inline constexpr uint64_t compute_mask(uint64_t const size) {
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << (size)) - 1u;
}


UnitigsDictionaryHash3::UnitigsDictionaryHash3() {}

UnitigsDictionaryHash3::UnitigsDictionaryHash3(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
}


inline size_t UnitigsDictionaryHash3::get_endpoints(const std::vector<std::vector<seqan3::dna4>> &input)
{
    size_t N = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        N += record.size();
        no_sequences++;
    }

    sdsl::bit_vector sequences = sdsl::bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = sdsl::sd_vector<>(sequences);
    endpoints_rank = sdsl::rank_support_sd<>(&endpoints);
    endpoints_select = sdsl::select_support_sd<>(&endpoints);

    return N;
} 


template <typename ViewType>
inline void UnitigsDictionaryHash3::mark_minimisers(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view, sdsl::bit_vector &r)
{
    for(auto & record : input) {
        for(auto && minimiser : record | view)
            r[minimiser.minimiser_value] = 1;
    }
}

template <typename ViewType>
inline void UnitigsDictionaryHash3::mark_unfreq_minimisers(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view, sdsl::bit_vector &r,
    uint8_t* count, sdsl::rank_support_v<1> &r_rank, uint8_t const m_thres)
{
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            if(count[i] < m_thres)
                r[minimiser.minimiser_value] = 1;
        }
    }
}

template <typename ViewType>
inline uint8_t* UnitigsDictionaryHash3::count_minimisers(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view,
    sdsl::rank_support_v<1> &r_rank, const uint8_t m_thres)
{
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    auto update_count = [&](uint64_t minimiser_value, size_t occurrences) {
        size_t i = r_rank(minimiser_value);
        count[i] += occurrences/span + 1;
        if(count[i] > m_thres)
            count[i] = m_thres;
    };

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view)
            update_count(minimiser.minimiser_value, minimiser.occurrences);
    }

    return count;
}


template <typename ViewType>
inline uint8_t* UnitigsDictionaryHash3::update_counts(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view, uint8_t* count,
    sdsl::bit_vector &r, sdsl::bit_vector &rnew, sdsl::rank_support_v<1> &r_rank, sdsl::rank_support_v<1> &rnew_rank)
{
    size_t c = rnew_rank(M);
    uint8_t* count_new = new uint8_t[c];
    std::memset(count_new, 0, c*sizeof(uint8_t));

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view)
            count_new[rnew_rank(minimiser.minimiser_value)] = count[r_rank(minimiser.minimiser_value)];
    }

    return count_new;
}

// template <int level>
// inline uint64_t UnitigsDictionaryHash3::build_S(uint8_t* count, sdsl::rank_support_v<1> &r_rank)
// {
//     size_t c = r_rank(M);
//     uint64_t n = 0;
//     for(uint64_t i=0; i < c; i++)
//         n += count[i];

//     sdsl::bit_vector s = sdsl::bit_vector(n+1, 0);
//     s[0] = 1;
//     size_t j = 0;
//     for (uint64_t i=0; i < c; i++) {
//         j += count[i];
//         s[j] = 1;
//     }

//     if constexpr (level == 1) {
//         s1 = s;
//         s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), n+1, 3);
//     }
//     if constexpr (level == 2) {
//         s2 = s;
//         s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), n+1, 3);
//     }

//     return n;
// }

inline uint64_t UnitigsDictionaryHash3::build_S1(uint8_t* count, sdsl::rank_support_v<1> &r_rank)
{
    size_t c = r_rank(M);
    uint64_t n = 0;
    for(uint64_t i=0; i < c; i++)
        n += count[i];

    s1 = sdsl::bit_vector(n+1, 0);
    s1[0] = 1;
    size_t j = 0;
    for (uint64_t i=0; i < c; i++) {
        j += count[i];
        s1[j] = 1;
    }

    s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), n+1, 3);

    return n;
}

template <typename ViewType>
inline void UnitigsDictionaryHash3::fill_offsets(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view,
    uint8_t* count, sdsl::bit_vector &r, sdsl::rank_support_v<1> &r_rank,
    sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> &s_select, int_vector<0> &offsets,
    uint64_t n, uint64_t N)
{
    offsets.width(std::bit_width(N));
    offsets.resize(n);

    size_t c = r_rank(M);
    std::memset(count, 0, c*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            if(r[minimiser.minimiser_value]) {
                size_t i = r_rank(minimiser.minimiser_value);
                size_t s = s_select.select(i);

                size_t o = minimiser.occurrences;
                uint64_t j = 0;
                while(o > span) {
                    offsets[s + count[i]] = length + minimiser.range_position + j*span;
                    count[i]++;
                    o -= span;
                    j++;
                }
                offsets[s + count[i]] = length + minimiser.range_position + j*span;
                count[i]++;
            }
        }
        length += sequence.size();
    }
}


// template <typename ViewType>
// template <int level>
// inline void UnitigsDictionaryHash3::compute_level(
//     const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view, const uint64_t text_length)
// {
//     if constexpr (level == 1) {
//         const uint8_t m_thres = m_thres1;
//         r1 = sdsl::bit_vector(M, 0);
//         sdsl::bit_vector r = r1;
//         sdsl::rank_support_v<1> r_rank = r1_rank;
//         sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select = s1_select;
//         int_vector<0> offsets = offsets1;
//     }
//     if constexpr (level == 2) {
//         const uint8_t m_thres = m_thres2;
//         r2 = sdsl::bit_vector(M, 0);
//         sdsl::bit_vector r = r2;
//         sdsl::rank_support_v<1> r_rank = r2_rank;
//         sux::bits::SimpleSelect<sux::util::AllocType::MALLOC> s_select = s2_select;
//         int_vector<0> offsets = offsets2;
//     }

//     std::cout << "filling bitvector R_tmp...\n";
//     sdsl::bit_vector rtmp = sdsl::bit_vector(M, 0);
//     mark_minimisers(input, view, rtmp);

//     std::cout << "count minimisers...\n";
//     sdsl::rank_support_v<1> rtmp_rank = sdsl::rank_support_v<1>(&rtmp);
//     uint8_t* count_tmp = count_minimisers(input, view, rtmp_rank, m_thres);

//     std::cout << "build bitvector R and according count array...\n";
//     mark_unfreq_minimisers(input, view, r, count_tmp, rtmp_rank, m_thres);
//     r_rank = sdsl::rank_support_v<1>(&r);
//     uint8_t* count = update_counts(input, view, count_tmp, rtmp, r, rtmp_rank, r_rank);

//     delete[] count_tmp;
//     // delete rtmp/rtmp_rank

//     std::cout << "filling bitvector S...\n";
//     uint64_t n = fill_S(count, r_rank);
//     std::cout << "filling offsets O...\n";
//     fill_offsets(input, view, count, r_rank, s_select, offsets, n, text_length);

//     delete[] count;
// }

template <typename ViewType>
inline void UnitigsDictionaryHash3::compute_level1(
    const std::vector<std::vector<seqan3::dna4>> &input, const ViewType& view, const uint64_t text_length)
{
    std::cout << "filling bitvector R_tmp...\n";
    sdsl::bit_vector rtmp = sdsl::bit_vector(M, 0);
    mark_minimisers(input, view, rtmp);

    std::cout << "count minimisers...\n";
    sdsl::rank_support_v<1> rtmp_rank = sdsl::rank_support_v<1>(&rtmp);
    uint8_t* count_tmp = count_minimisers(input, view, rtmp_rank, m_thres1);

    std::cout << "build bitvector R...\n";
    sdsl::bit_vector r1 = sdsl::bit_vector(M, 0);
    mark_unfreq_minimisers(input, view, r1, count_tmp, rtmp_rank, m_thres1);
    r1_rank = sdsl::rank_support_v<1>(&r1);
    std::cout << "build according count array...\n";
    uint8_t* count = update_counts(input, view, count_tmp, rtmp, r1, rtmp_rank, r1_rank);

    delete[] count_tmp;
    // delete rtmp/rtmp_rank

    std::cout << "filling bitvector S...\n";
    // uint64_t n = build_S<1>(count, r1_rank);
    uint64_t n = build_S1(count, r1_rank);
    std::cout << "filling offsets O...\n";
    fill_offsets(input, view, count, r1, r1_rank, s1_select, offsets1, n, text_length);

    delete[] count;
}


template <typename ViewType>
inline std::vector<std::vector<seqan3::dna4>> UnitigsDictionaryHash3::get_frequent_skmers(
    const std::vector<std::vector<seqan3::dna4>> &input,
    const ViewType& view, sdsl::bit_vector &r)
{
    std::vector<std::vector<seqan3::dna4>> freq_skmers;
    for(auto & sequence : input) {
        size_t start_position;
        bool freq;
        bool current_freq;
        for(auto && minimiser : sequence | view) {
            freq = !r[minimiser.minimiser_value];
            if(freq)
                start_position = 0;
            break;
        }
        for(auto && minimiser : sequence | view) {
            current_freq = !r[minimiser.minimiser_value];

            if(!freq && current_freq)
                start_position = minimiser.range_position;
            if(freq && !current_freq) {
                std::vector<seqan3::dna4> freq_skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    freq_skmer.push_back(sequence[i]);
                freq_skmers.push_back(freq_skmer);
            }

            freq = current_freq;
        }
    }

    return freq_skmers;
}


int UnitigsDictionaryHash3::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    M = 1ULL << (m+m);
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});

    std::cout << "get endpoints...\n";
    uint64_t text_length = get_endpoints(input);

    std::cout << "compute level 1...\n";
    // compute_level<1>(input, view1, text_length);
    compute_level1(input, view1, text_length);

    std::cout << "extract uncovered sequence parts...\n";
    std::vector<std::vector<seqan3::dna4>> remaining_skmers1;
    remaining_skmers1 = get_frequent_skmers(input, view1, r1);

    size_t len_rem_seqs = 0;
    for(auto & sequence : remaining_skmers1)
        len_rem_seqs += sequence.size();

    // std::cout << "remaining superkmers " << remaining_skmers1.size() << " (" << (double)remaining_skmers1.size()/skmers*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/text_length*100 << "%)\n";

    // auto view2 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed2});

    // std::cout << "compute level 2...\n";
    // compute_level(remaining_skmers1, view2, 2);

    // std::cout << "extract uncovered sequence parts...\n";
    // std::vector<std::vector<seqan3::dna4>> remaining_sequences2;
    // remaining_sequences2 = get_frequent_skmers(remaining_sequences1, view2, r2);

    // len_rem_seqs = 0;
    // for(auto & sequence : remaining_sequences2)
    //     len_rem_seqs += sequence.size();
    // std::cout << "remaining superkmers " << remaining_sequences2.size() << " (" << (double)remaining_sequences2.size()/skmers*100 << "%) ";
    // std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";

    // std::cout << "filling HT...\n";
    // std::unordered_set<uint64_t> freq_kmers;
    // std::unordered_set<uint64_t> freq_minimzer;

    // auto view3 = srindex::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k, .seed=seed2});
    // for(auto & sequence : remaining_sequences2) {
    //     for(auto && minimiser : sequence | view3) {
    //         freq_minimzer.insert(minimiser.minimiser_value);
    //         // freq_kmers.insert(minimiser.window_value);
    //     }
    // }

    std::cout << "copy text...\n";
    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

    // std::cout << "====== report ======\n";
    // std::cout << "text length: " << text_length << "\n";
    // std::cout << "no kmers: " << kmers <<  '\n';
    
    // std::cout << "no minimiser: " << n << "\n";
    // std::cout << "no distinct minimiser: " << c << "\n";
    // std::cout << "minimiser going to level 2: " << c-c1 << "  " << (double) (c-c1)/c*100 << "%\n";
    // std::cout << "no minimiser1: " << n1 << "  " << (double) n1/(n1+n2)*100 << "%\n";
    // std::cout << "no distinct minimiser1: " << c1 << "  " << (double) c1/(c1+c2)*100 << "%\n";
    // std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    // std::cout << "no minimiser2: " << n2 << "  " << (double) n2/(n1+n2)*100 << "%\n";
    // std::cout << "no distinct minimiser2: " << c2 << "  " << (double) c2/(c1+c2)*100 << "%\n";
    // std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';

    // std::cout << "density r1: " << (double) c1/M*100 << "%\n";
    // std::cout << "density r2: " << (double) c2/M*100 << "%\n";
    // std::cout << "density s1: " << (double) s1_select.bitCount()/(n1+1)*100 <<  "%\n";
    // std::cout << "density s2: " << (double) s2_select.bitCount()/(n2+1)*100 <<  "%\n";
    // std::cout << "\nspace per kmer in bit:\n";
    // std::cout << "text: " << (double) 2*N/kmers << "\n";
    // std::cout << "endpoints: " << (double) size_in_bytes(endpoints)/(8*kmers) << "\n";
    // std::cout << "offsets1: " << (double) n1*offset_width/kmers << "\n";
    // std::cout << "offsets2: " << (double) n2*offset_width/kmers << "\n";
    // std::cout << "R_1: " << (double) M/kmers << "\n";
    // std::cout << "R_2: " << (double) M/kmers << "\n";
    // std::cout << "S_1: " << (double) (n1+1)/kmers << "\n";
    // std::cout << "S_2: " << (double) (n2+1)/kmers << "\n";

    // std::cout << "total: " << (double) (n1*offset_width+n2*offset_width+2*N+M+M+n1+1+n2+1+size_in_bytes(endpoints)/8)/kmers << "\n";

    return 0;
}

template <int level>
inline void UnitigsDictionaryHash3::fill_buffer(std::vector<uint64_t> &buffer, const uint64_t mask, size_t p, size_t q)
{
    for(uint64_t i = 0; i < q-p; i++) {
        uint64_t hash = 0;
        size_t o;
        if constexpr (level == 1)
            o = offsets1[p+i];
        if constexpr (level == 2)
            o = offsets2[p+i];

        size_t next_endpoint = endpoints_select(endpoints_rank(o+1)+1);
        size_t e = o+k+span;
        if(e > next_endpoint)
            e = next_endpoint;
        for (uint64_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
        }
        buffer.push_back(hash);
        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buffer.push_back(hash);
        }
    }
}


// inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query, uint64_t query_rev) {
//     for(uint64_t i=0; i < array.size(); i++) {
//         if(array[i] == query || array[i] == query_rev)
//             return true;
//     }
//     return false;
// }

inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, uint64_t &last_found) {
    for(uint64_t i=last_found+1; i < array.size(); i++) {
        if(array[i] == query || array[i] == queryrc) {
            last_found = i;
            return true;
        }
    }
    for(uint64_t i=0; i < last_found+1; i++) {
        if(array[i] == query || array[i] == queryrc) {
            last_found = i;
            return true;
        }
    }
    return false;
}


uint64_t UnitigsDictionaryHash3::streaming_query(const std::vector<seqan3::dna4> &query)
{
    auto view = srindex::views::two_minimisers_and_window_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});

    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(2u * k);
    uint64_t current_minimiser = 0; // lets hope the first minimiser is not 0
    std::vector<uint64_t> buffer;
    uint64_t last_found = 0;

    for(auto && minimisers : query | view)
    {
        if(minimisers.minimiser1_value == current_minimiser) {
            occurences += lookup_serial(buffer, minimisers.window_value, minimisers.window_value_rev, last_found);
        }
        else {
            if(r1[minimisers.minimiser1_value]) {
                size_t minimizer_id = r1_rank(minimisers.minimiser1_value);
                size_t p = s1_select.select(minimizer_id);
                size_t q = s1_select.select(minimizer_id+1);

                buffer.clear();
                last_found = 0;
                // fill_buffer1(buffer, mask, p, q);
                fill_buffer<1>(buffer, mask, p, q);
                occurences += lookup_serial(buffer, minimisers.window_value, minimisers.window_value_rev, last_found);
                current_minimiser = minimisers.minimiser1_value;
            }
            // else {
            //     if(r2[minimisers.minimiser1_value]) {
            //         size_t minimizer_id = r2_rank(minimisers.minimiser1_value);
            //         size_t p = s2_select.select(minimizer_id);
            //         size_t q = s2_select.select(minimizer_id+1);

            //         buffer.clear();
            //         last_found = 0;
            //         // fill_buffer2(buffer, mask, p, q);
            //         fill_buffer<2>(buffer, mask, p, q);
            //         occurences += lookup_serial(buffer, minimisers.window_value, minimisers.window_value_rev, last_found);
            //         current_minimiser = minimisers.minimiser1_value;
            //     }
            // }
        }
    }

    return occurences;
}


int UnitigsDictionaryHash3::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(r1, out);
    // seqan3::contrib::sdsl::serialize(r2, out);
    seqan3::contrib::sdsl::serialize(s1, out);
    // seqan3::contrib::sdsl::serialize(s2, out);
    seqan3::contrib::sdsl::serialize(this->offsets1, out);
    // seqan3::contrib::sdsl::serialize(this->offsets2, out);
    seqan3::contrib::sdsl::serialize(this->endpoints, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->text);

    out.close();
    return 0;
}

int UnitigsDictionaryHash3::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::load(r1, in);
    r1_rank = rank_support_v<1>(&r1);
    // seqan3::contrib::sdsl::load(r2, in);
    // r2_rank = rank_support_v<1>(&r2);
    seqan3::contrib::sdsl::load(s1, in);
    this->s1_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s1.data()), s1.size(), 3);
    // seqan3::contrib::sdsl::load(s2, in);
    // this->s2_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s2.data()), s2.size(), 3);
    seqan3::contrib::sdsl::load(this->offsets1, in);
    // seqan3::contrib::sdsl::load(this->offsets2, in);
    seqan3::contrib::sdsl::load(this->endpoints, in);
    endpoints_rank = rank_support_sd<>(&endpoints);
    endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);

    in.close();
    return 0;
}

