#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>

#include "dict.h"
#include "kmer_minimiser_hash.hpp"



LookupUnitigsDictionary::LookupUnitigsDictionary() {}

LookupUnitigsDictionary::LookupUnitigsDictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
}


int LookupUnitigsDictionary::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    auto view = bsc::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    const uint64_t M = 1 << (m+m); // 4^m

    // todo: simons bitvector
    r = bit_vector(M);
    // todo: faster minimiser only view
    size_t N = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r[minimiser.minimiser_value] = 1;
        }
        N += record.size();
    }
    r_rank = rank_support_v<1>(&r);

    size_t c = r_rank(M);
    std::cout << "no distinct kmers: " << c << '\n';
    std::cout << "density r: " << (float) c/M*100 << "%\n";

    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    // todo: init hashtable with adequate size, faster hashtable
    std::unordered_map<uint64_t, uint32_t> cb;
    
    uint64_t n = 0;
    uint64_t kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            uint64_t i = r_rank(minimiser.minimiser_value);
            uint64_t o = minimiser.occurrences;

            int w = o/k + 1;
            if(count[i] == 255)
                cb[i] += w;
            else if(count[i] + w >= 255) {
                cb[i] = w - (255 - count[i]);
                count[i] = 255;
            }
            else
                count[i] += w;
            n += w;
            kmers += o;
        }
    }
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "freq minimiser: " << cb.size() << '\n';

    s = bit_vector(n+1, 0);
    s[0] = 1;
    int j = 0;
    for (size_t i=0; i < c; i++) {
        j += count[i];
        if(count[i] == 255)
            j += cb[i];
        s[j] = 1;
    }
    simple_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), n+1, 3);

    std::cout << "density s: " << (float) simple_select.bitCount()/(n+1)*100 <<  "%\n";

    offset_width = std::bit_width(N);
    span_width = std::bit_width(k);
    offset.width(offset_width);
    offset.resize(n);
    span.width(span_width);
    span.resize(n);

    std::cout << "offset width: " << offset_width << "\n";
    std::cout << "span width: " << span_width << "\n";
    std::cout << "allocated " << n*offset_width/8 << " bytes for O\n";

    std::memset(count, 0, c*sizeof(uint8_t));
    // faster?
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    uint64_t cur_text_length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view) {
            uint64_t i = r_rank(minimiser.minimiser_value);
            uint64_t s = simple_select.select(i);
            int o = minimiser.occurrences;
            int j = 0;
            while(o > k) {
                if(count[i] == 255) {
                    offset[s + count[i] + cb[i]] = cur_text_length + minimiser.range_position + j*k;
                    span[s + count[i] + cb[i]] = k;
                    cb[i]++;
                }
                else {
                    offset[s + count[i]] = cur_text_length + minimiser.range_position + j*k;
                    span[s + count[i]] = k;
                    count[i]++;
                }
                o -= k;
                j++;
            }
            if(count[i] == 255) {
                offset[s + count[i] + cb[i]] = cur_text_length + minimiser.range_position + j*k;
                span[s + count[i] + cb[i]] = o;
                cb[i]++;
            }
            else {
                offset[s + count[i]] = cur_text_length + minimiser.range_position + j*k;
                span[s + count[i]] = o;
                count[i]++;
            }
        }
        cur_text_length += sequence.size();
    }

    delete[] count; // delete cb?


    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

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


int LookupUnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &query)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

    int occurences = 0;
    const uint64_t mask = compute_mask(k);

    for(auto && minimiser : query | query_view) {
        if(r[minimiser.minimiser_value]) {
            size_t r = r_rank(minimiser.minimiser_value);
            size_t p = simple_select.select(r);
            size_t q = simple_select.select(r+1);
            size_t b = q - p;

            bool found = false;
            for(int i = 0; !found && i < b; i++) {
                uint64_t hash = 0;
                size_t o = offset[p+i];
                for (int j=o; j < o+k; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                }
                if(minimiser.window_value == hash) {
                    occurences++;
                    break;
                }
                for(int j=o+k; j < o+k+span[p+i]; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                    if(minimiser.window_value == hash) {
                        found = true;
                        occurences++;
                        break;
                    }
                }
            }
        }
        // todo: else ht lookup
    }

    return occurences;
}


// int LookupUnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &text,
//     const std::vector<seqan3::dna4> &query, std::vector<std::pair<uint64_t, uint64_t>> &positions)
// {
//     auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

//     const size_t N = text.size();
//     const uint64_t mask = compute_mask(k);

//     int kmer = 0;
//     for(auto && minimiser : query | query_view) {
//         if(r[minimiser.minimiser_value]) {
//             size_t i = r_rank(minimiser.minimiser_value);
//             size_t p = simple_select.select(i);
//             size_t q = simple_select.select(i+1);
//             size_t b = q - p;

//             bool found = false;
//             for(int i = 0; !found && i < b; i++) {
//                 size_t o = offset[p+i];
//                 size_t e = std::min(o+2*k, N);

//                 uint64_t hash = 0;
//                 for (int j=o; j < o+k; j++) {
//                     uint64_t const new_rank = seqan3::to_rank(text[j]);
//                     hash <<= 2;
//                     hash |= new_rank;
//                     hash &= mask;
//                 }
//                 if(minimiser.window_value == hash) {
//                     positions.push_back({kmer, o});
//                     break;
//                 }
//                 for(int j=o+k; j < e; j++) {
//                     uint64_t const new_rank = seqan3::to_rank(text[j]);
//                     hash <<= 2;
//                     hash |= new_rank;
//                     hash &= mask;
//                     if(minimiser.window_value == hash) {
//                         found = true;
//                         positions.push_back({kmer, j-k+1});
//                         break;
//                     }
//                 }
//             }
//         }
//         kmer++;
//     }

//     return 0;
// }


int LookupUnitigsDictionary::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(r, out); // save r_rank?
    seqan3::contrib::sdsl::serialize(s, out);
    seqan3::contrib::sdsl::serialize(this->offset, out);
    seqan3::contrib::sdsl::serialize(this->span, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->text);

    out.close();
    return 0;
}

int LookupUnitigsDictionary::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::load(r, in);
    r_rank = rank_support_v<1>(&r);
    seqan3::contrib::sdsl::load(s, in);
    this->simple_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), s.size(), 3);
    seqan3::contrib::sdsl::load(this->offset, in);
    seqan3::contrib::sdsl::load(this->span, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);

    in.close();
    return 0;
}

