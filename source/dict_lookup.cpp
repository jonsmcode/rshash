#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "dict.h"
#include "kmer_minimiser_hash.hpp"


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


LookupDictionary::LookupDictionary() {}

LookupDictionary::LookupDictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
}


int LookupDictionary::build(const std::vector<seqan3::dna4> &text)
{
    auto view = bsc::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    const uint64_t M = 1 << (m+m); // 4^m

    r = bit_vector(M);
    for(auto && minimiser : text | view) {
        r[minimiser.minimiser_value] = 1;
    }
    r_rank = rank_support_v<1>(&r);

    uint64_t c = r_rank(M);
    std::cout  << c << "\n";

    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    // todo: init hashtable with adequate size, faster hashtable
    std::unordered_map<uint64_t, uint32_t> cb;
    
    int n = 0;
    // (todo: implement minmizer s.t. minmizer range at most window length)
    for(auto && minimiser : text | view) {
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
    }
    std::cout << "no minimiser: " << n <<  "\n";
    std::cout << "freq kmers: " << (float) cb.size()/n *100 << "%\n";

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

    size_t offset_width = std::bit_width(text.size());
    offset.width(offset_width);
    offset.resize(n);

    std::cout << "allocated memory.\n";

    std::memset(count, 0, c*sizeof(uint8_t));
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    for (auto && minimiser : text | view) {
        uint64_t i = r_rank(minimiser.minimiser_value);
        uint64_t s = simple_select.select(i);
        int o = minimiser.occurrences;
        int j = 0;
        while(o > k) {
            if(count[i] == 255) {
                offset[s + count[i] + cb[i]] = minimiser.range_position + j*k;
                cb[i]++;
            }
            else {
                offset[s + count[i]] = minimiser.range_position + j*k;
                count[i]++;
            }
            o -= k;
            j++;
        }
        if(count[i] == 255) {
            offset[s + count[i] + cb[i]] = minimiser.range_position + j*k;
            cb[i]++;
        }
        else {
            offset[s + count[i]] = minimiser.range_position + j*k;
            count[i]++;
        }
    }

    delete[] count; // delete cb?

    return 0;
}

int LookupDictionary::streaming_query(const std::vector<seqan3::dna4> &text,
                                      const std::vector<seqan3::dna4> &query)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});
    uint64_t const seed = 0;
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
                      | std::views::transform([](uint64_t i) {return i ^ seed;});

    int occurences = 0;

    for(auto && minimiser : query | query_view) {
        if(r[minimiser.minimiser_value]) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t p = simple_select.select(i);
            size_t q = simple_select.select(i+1);
            size_t b = q - p;

            bool found = false;
            for(int i = 0; !found && i < b; i++) {
                for (auto && hash : text | std::views::drop(offset[p+i])
                                         | std::views::take(k - m + 1) | kmer_view) {
                    if(minimiser.window_value == hash) {
                        found = true;
                        occurences++;
                        break;
                    }
                }
            }
        }
    }

    return occurences;
}

int LookupDictionary::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::rrr_vector<> rrr(r);
    seqan3::contrib::sdsl::serialize(rrr, out);
    seqan3::contrib::sdsl::sd_vector<> sds(s);
    seqan3::contrib::sdsl::serialize(sds, out);
    seqan3::contrib::sdsl::serialize(this->offset, out);
    out.close();
    return 0;
}

int LookupDictionary::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::rrr_vector<> rrr;
    seqan3::contrib::sdsl::load(rrr, in);
    this->r = bit_vector(rrr.size());
    for (size_t i = 0; i < rrr.size(); i++)
        r[i] = rrr[i];
    r_rank = rank_support_v<1>(&r);
    seqan3::contrib::sdsl::sd_vector<> sds;
    seqan3::contrib::sdsl::load(sds, in);
    this->s = bit_vector(sds.size());
    for (size_t i = 0; i < s.size(); i++)
        s[i] = sds[i];
    this->simple_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), s.size(), 3);
    seqan3::contrib::sdsl::load(this->offset, in);
    in.close();
    return 0;
}

