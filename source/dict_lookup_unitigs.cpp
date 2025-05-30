#include <filesystem>
// #include <cstdint>
// #include <immintrin.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>

#include <cereal/archives/binary.hpp>

#include "dict.h"
#include "kmer_minimiser_hash.hpp"


const uint8_t m_thres = 5;



UnitigsDictionary::UnitigsDictionary() {}

UnitigsDictionary::UnitigsDictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
}


int UnitigsDictionary::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view = bsc::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    const uint64_t M = 1 << (m+m); // 4^m

    std::cout << "extracting minimizers...\n";
    r = bit_vector(M, 0);

    size_t N = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r[minimiser.minimiser_value] = 1;
        }
        N += record.size();
        no_sequences++;
    }
    r_rank = rank_support_v<1>(&r);

    bit_vector sequences = bit_vector(N+1, 0);
    uint64_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = seqan3::contrib::sdsl::sd_vector<>(sequences);
    endpoints_rank = seqan3::contrib::sdsl::rank_support_sd<>(&endpoints);
    endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);

    std::cout << "counting minimizers...\n";
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));
    std::unordered_map<uint64_t, uint32_t> cb;
    uint32_t max_minimizer_occs = 0;

    uint64_t n = 0;
    uint64_t kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            int w = o/k + 1;
            if(count[i] == m_thres) {
                cb[i] += w;
                max_minimizer_occs = std::max(cb[i], max_minimizer_occs);
            }
            else if(count[i] + w >= m_thres) {
                cb[i] = w - (m_thres - count[i]);
                count[i] = m_thres;
                max_minimizer_occs = std::max(cb[i], max_minimizer_occs);
            }
            else {
                count[i] += w;
                // max_minimizer_occs = std::max(count[i], max_minimizer_occs);
            }
            n += w;
            kmers += o;
        }
    }

    std::cout << "max minimizer occurences: " << max_minimizer_occs << "\n";

    uint64_t freq_kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            if(count[i] == m_thres)
                freq_kmers++;
        }
    }
    std::cout << "freq k-mers: " << freq_kmers << "\n";
    

    // for (const auto& [m, occ] : cb) {
    //     r[m] = 0;
    // }

    std::cout << "filling bitvector S...\n";
    s = bit_vector(n+1, 0);
    s[0] = 1;
    j = 0;
    for (uint64_t i=0; i < c; i++) {
        j += count[i];
        if(count[i] == m_thres)
            j += cb[i];
        s[j] = 1;
    }
    s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), n+1, 3);

    size_t offset_width = std::bit_width(N);
    offsets.width(offset_width);
    offsets.resize(n);
    std::cout << "allocated " << n*(offset_width)/8000000 << " Mb for O\n";

    std::memset(count, 0, c*sizeof(uint8_t));
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    std::cout << "filling O...\n";
    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t s = s_select.select(i);
            size_t o = minimiser.occurrences;
            uint64_t j = 0;
            while(o > k) {
                if(count[i] == m_thres) {
                    offsets[s + count[i] + cb[i]] = length + minimiser.range_position + j*k;
                    cb[i]++;
                }
                else {
                    offsets[s + count[i]] = length + minimiser.range_position + j*k;
                    count[i]++;
                }
                o -= k;
                j++;
            }
            if(count[i] == m_thres) {
                offsets[s + count[i] + cb[i]] = length + minimiser.range_position + j*k;
                cb[i]++;
            }
            else {
                offsets[s + count[i]] = length + minimiser.range_position + j*k;
                count[i]++;
            }
        }
        length += sequence.size();
    }

    delete[] count;

    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no distinct kmers: " << c << '\n';
    std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "freq minimiser: " << cb.size() << ", " << (float) cb.size()/n*100 << "%\n";
    std::cout << "density r: " << (double) c/M*100 << "%\n";
    std::cout << "density s: " << (double) s_select.bitCount()/(n+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    // std::cout << "endpoints: " << (double) std::bit_width(max_seq_length)*no_sequences/(8*kmers) << "\n";
    // std::cout << "endpoints: " << (double) endpoints.size()/(8*kmers) << "\n";
    std::cout << "endpoints: " << (double) size_in_bytes(endpoints)/(8*kmers) << "\n";
    std::cout << "offsets: " << (double) n*offset_width/kmers << "\n";
    std::cout << "R: " << (double) M/kmers << "\n";
    std::cout << "S: " << (double) (n+1)/kmers << "\n";
    std::cout << "HT: " << (double) sizeof(*cbk.begin())*cbk.size()/(8*kmers) << "\n";

    std::cout << "total: " << (double) (n*offset_width+2*N+M+n+1+size_in_bytes(endpoints)/8)/kmers << "\n";

    return 0;
}


static inline constexpr uint64_t compute_mask(uint64_t const kmer_size)
{
    assert(kmer_size > 0u);
    assert(kmer_size <= 32u);

    if (kmer_size == 32u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << (2u * kmer_size)) - 1u;
}


inline void UnitigsDictionary::fill_buffer(std::vector<uint64_t> &buffer, const uint64_t mask, size_t p, size_t q) {
    for(uint64_t i = 0; i < q-p; i++) {
        uint64_t hash = 0;
        size_t o = offsets[p+i];
        for (uint64_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
        }
        buffer.push_back(hash);
        size_t next_endpoint = endpoints_select(endpoints_rank(o+1)+1);
        size_t e = o+k+k;
        if(e > next_endpoint)
            e = next_endpoint;
        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buffer.push_back(hash);
        }
    }
}


inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query) {
    for(uint64_t i=0; i < array.size(); i++) {
        if(array[i] == query)
            return true;
    }
    return false;
}

// inline bool contains_uint64_avx512(const std::vector<uint64_t>& arr, uint64_t value) {
//     const uint64_t* data = arr.data();
//     size_t len = arr.size();
//     __m512i target = _mm512_set1_epi64(value);
//     size_t i = 0;
//     for (; i + 8 <= len; i += 8) {
//         __m512i chunk = _mm512_loadu_si512((__m512i const*)(data + i));
//         __mmask8 mask = _mm512_cmpeq_epi64_mask(chunk, target);
//         if (mask != 0)
//             return true;
//     }
//     for (; i < len; ++i) {
//         if (data[i] == value)
//             return true;
//     }
//     return false;
// }


uint64_t UnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &query)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(k);
    uint64_t current_minimiser = 0; // lets hope the first minimiser is not 0
    std::vector<uint64_t> buffer;

    for(auto && minimiser : query | query_view)
    {
        if(minimiser.minimiser_value == current_minimiser) {
            occurences += lookup_serial(buffer, minimiser.window_value);
        }
        else {
            if(r[minimiser.minimiser_value]) {
                size_t minimizer_id = r_rank(minimiser.minimiser_value);
                size_t p = s_select.select(minimizer_id);
                size_t q = s_select.select(minimizer_id+1);

                buffer.clear();
                fill_buffer(buffer, mask, p, q);
                occurences += lookup_serial(buffer, minimiser.window_value);
                current_minimiser = minimiser.minimiser_value;
            }
            // else {
            //    occurences += cbk.contains(minimiser.window_value);
            // }
        }
        
    }

    return occurences;
}



inline void UnitigsDictionary::fill_buffer(std::vector<std::vector<uint64_t>> &buffer, std::vector<std::vector<uint64_t>> &positions,
                                           std::vector<uint64_t> &sequences, const uint64_t mask, size_t p, size_t q)
{
    for(uint64_t i = 0; i < q-p; i++)
    {
        std::vector<uint64_t> buff;
        std::vector<uint64_t> pos;

        uint64_t hash = 0;
        size_t o = offsets[p+i];
        size_t sequence_id = endpoints_rank(o+1);
        size_t endpoint = endpoints_select(sequence_id);
        size_t next_endpoint = endpoints_select(sequence_id+1);

        size_t e = o+k+k;
        if(e > next_endpoint)
            e = next_endpoint;

        for (uint64_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
        }
        buff.push_back(hash);
        pos.push_back(o-endpoint);
        for(size_t j=o+k; j < e; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
            hash &= mask;
            buff.push_back(hash);
            pos.push_back(j-k+1-endpoint);
        }

        buffer.push_back(buff);
        positions.push_back(pos);
        sequences.push_back(sequence_id);
    }
}


inline void locate_serial(std::vector<std::vector<uint64_t>> &buffer, std::vector<std::vector<uint64_t>> &positions,
                          std::vector<uint64_t> &sequences, const uint64_t query, const uint64_t kmer,
                          std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &result) {
    for(uint64_t i=0; i < buffer.size(); i++) {
        for(uint64_t j=0; j < buffer[i].size(); j++) {
            if(buffer[i][j] == query)
                result.push_back({kmer, sequences[i], positions[i][j]});
        }
    }
}


uint64_t UnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &query,
                                            std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &result)
{
    auto view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

    const uint64_t mask = compute_mask(k);
    uint64_t kmer = 0;
    uint64_t current_minimiser = 0;
    std::vector<std::vector<uint64_t>> buffer;
    std::vector<std::vector<uint64_t>> positions;
    std::vector<uint64_t> sequences;

    for(auto && minimiser : query | view) {
        if(minimiser.minimiser_value == current_minimiser) {
            locate_serial(buffer, positions, sequences, minimiser.window_value, kmer, result);
        }
        else {
            if(r[minimiser.minimiser_value]) {
                size_t minimizer_id = r_rank(minimiser.minimiser_value);
                size_t p = s_select.select(minimizer_id);
                size_t q = s_select.select(minimizer_id+1);

                buffer.clear();
                positions.clear();
                sequences.clear();
                fill_buffer(buffer, positions, sequences, mask, p, q);
                locate_serial(buffer, positions, sequences, minimiser.window_value, kmer, result);

                current_minimiser = minimiser.minimiser_value;
            }
        }
        kmer++;

    }

    return 0;
}


int UnitigsDictionary::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(r, out); // save r_rank
    seqan3::contrib::sdsl::serialize(s, out);
    // out << this->s_select;
    seqan3::contrib::sdsl::serialize(this->offsets, out);
    seqan3::contrib::sdsl::serialize(this->endpoints, out);


    cereal::BinaryOutputArchive archive(out);
    archive(this->text);

    out.close();
    return 0;
}

int UnitigsDictionary::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::load(r, in);
    r_rank = rank_support_v<1>(&r);
    seqan3::contrib::sdsl::load(s, in);
    // this->s_select << in;
    this->s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), s.size(), 3);
    seqan3::contrib::sdsl::load(this->offsets, in);
    seqan3::contrib::sdsl::load(this->endpoints, in);
    endpoints_rank = rank_support_sd<>(&endpoints);
    endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);

    in.close();
    return 0;
}

