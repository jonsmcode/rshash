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


Dictionary::Dictionary() {}

Dictionary::Dictionary(uint8_t const k, uint8_t const m) {
    this->k = k;
    this->m = m;
}


int Dictionary::build(const std::vector<seqan3::dna4> &text)
{
    auto view = bsc::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    const uint64_t M = 1 << (m+m); // 4^m

    r = bit_vector(M);
    for(auto && minimiser : text | view) {
        r[minimiser.minimiser_value] = 1;
    }
    // todo: dense bv, faster library
    sdr = sd_vector<>(r);
    r_rank = sd_vector<>::rank_1_type(&sdr);

    uint64_t c = r_rank(M);
    std::cout  << c <<  "\n";

    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    // todo: init hashtable with adequate size, faster hashtable
    std::unordered_map<uint64_t, uint32_t> cb;
    
    int n = 0;
    // (todo: implement minmizer s.t. minmizer range at most window length)
    for(auto && minimiser : text | view) {
        uint64_t r = r_rank(minimiser.minimiser_value);
        uint64_t o = minimiser.occurrences;

        int w = o/k + 1;
        if(count[r] == 255)
            cb[r] += w;
        else if(count[r] + w >= 255) {
            cb[r] = w - (255 - count[r]);
            count[r] = 255;
        }
        else
            count[r] += w;
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

    // todo: test runtime allocated (and compressed) sdsl::int_vector<0> sizes
    sds = sd_vector<>(s);
    s_select = sd_vector<>::select_1_type(&sds);

    size_t offset_width = std::bit_width(text.size());
    size_t span_width = std::bit_width(k);
    // offset.width(offset_width + span_width);
    offset.width(offset_width);
    offset.resize(n);
    span.width(span_width);
    span.resize(n);

    std::cout << "allocated memory.\n";
    // std::cout << "offset width " << offset_width << ", span width " << span_width << std::endl;
    // std::cout << "total width " << (int) offset.width() << ", total bits " << offset.bit_size() << std::endl;

    std::memset(count, 0, c*sizeof(uint8_t));
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    for (auto && minimiser : text | view) {
        // seqan3::debug_stream << minimiser.minimiser_value << ',' << kmer_to_string(minimiser.minimiser_value, m)
                             // << ',' << minimiser.range_position << ',' << minimiser.occurrences << '\n';
        uint64_t r = r_rank(minimiser.minimiser_value);
        uint64_t s = s_select(r+1);
        // seqan3::debug_stream << r << ' ' << s << '\n';
        int o = minimiser.occurrences;
        int j = 0;
        while(o > k) {
            // offset[s + count[r]] = (minimiser.range_position + j*k << offset_width) | k;
            if(count[r] == 255) {
                offset[s + count[r] + cb[r]] = minimiser.range_position + j*k;
                span[s + count[r] + cb[r]] = k;
                cb[r]++;
            }
            else {
                offset[s + count[r]] = minimiser.range_position + j*k;
                span[s + count[r]] = k;
                count[r]++;
            }
            o -= k;
            j++;
        }
        // offset[s + count[r]] = (minimiser.range_position + j*k << offset_width) | o;
        if(count[r] == 255) {
            offset[s + count[r] + cb[r]] = minimiser.range_position + j*k;
            span[s + count[r] + cb[r]] = o;
            cb[r]++;
        }
        else {
            offset[s + count[r]] = minimiser.range_position + j*k;
            span[s + count[r]] = o;
            count[r]++;
        }
    }

    delete[] count; // delete cb?

    return 0;
}


int Dictionary::streaming_query(const std::vector<seqan3::dna4> &text,
                                const std::vector<seqan3::dna4> &query,
                                std::vector<uint64_t> &positions)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});
    // uint64_t const seed = 347692;
    uint64_t const seed = 0;
    // uint64_t const seed = 0x8F'3F'73'B5'CF'1C'9A'DE;
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
                      | std::views::transform([](uint64_t i) {return i ^ seed;});

    // todo: test hashtable instead of vector
    std::vector<uint64_t> kmerbuffer;
    std::vector<uint64_t> startpositions;
    int cur_minimiser = -1;

    for(auto && minimiser : query | query_view) {
        if(sdr[minimiser.minimiser_value]) {
            if(minimiser.minimiser_value != cur_minimiser) {
                kmerbuffer.clear();
                startpositions.clear();

                // todo: fast select bitvector implementation
                size_t r = r_rank(minimiser.minimiser_value);
                size_t p = s_select(r+1);
                size_t q = s_select(r+2);
                size_t b = q - p;

                for(int i = 0; i < b; i++) {
                    // todo: test with(out) span
                    // size_t sp = k - m + 1;
                    size_t sp = span[p+i]+k-1;

                    int j = 0;
                    for (auto && hash : text | std::views::drop(offset[p+i])
                                             | std::views::take(sp) | kmer_view) {
                        // seqan3::debug_stream << kmer_to_string(hash, k) << '\n';
                        kmerbuffer.push_back(hash);
                        startpositions.push_back(offset[p+i] + j);
                        j++;
                    }
                }

                cur_minimiser = minimiser.minimiser_value;
            }
            // todo: simd for buffer
            for(int i = 0; i < kmerbuffer.size(); i++) {
                if(minimiser.window_value == kmerbuffer[i])
                    positions.push_back(startpositions[i]);
            }
        }
        
    }

    return 0;
}


int Dictionary::streaming_query(const std::vector<seqan3::dna4> &text,
                                const std::vector<seqan3::dna4> &query)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});
    uint64_t const seed = 0;
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
                      | std::views::transform([](uint64_t i) {return i ^ seed;});

    std::vector<uint64_t> kmerbuffer;
    int cur_minimiser = -1;
    int occurences = 0;

    for(auto && minimiser : query | query_view) {
        if(sdr[minimiser.minimiser_value]) {
            if(minimiser.minimiser_value != cur_minimiser) {
                kmerbuffer.clear();

                size_t r = r_rank(minimiser.minimiser_value);
                size_t p = s_select(r+1);
                size_t q = s_select(r+2);
                size_t b = q - p;

                for(int i = 0; i < b; i++) {
                    size_t sp = span[p+i]+k-1;

                    int j = 0;
                    for (auto && hash : text | std::views::drop(offset[p+i])
                                             | std::views::take(sp) | kmer_view) {
                        kmerbuffer.push_back(hash);
                        j++;
                    }
                }
                cur_minimiser = minimiser.minimiser_value;
            }

            for(int i = 0; i < kmerbuffer.size(); i++) {
                if(minimiser.window_value == kmerbuffer[i]) {
                    occurences++;
                    break;
                }
            }
        }

    }

    return occurences;
}


int Dictionary::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(this->sdr, out);
    seqan3::contrib::sdsl::serialize(this->sds, out);
    seqan3::contrib::sdsl::serialize(this->offset, out);
    seqan3::contrib::sdsl::serialize(this->span, out);
    out.close();
    return 0;
}

int Dictionary::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::load(this->sdr, in);
    this->r_rank = sd_vector<>::rank_1_type(&sdr);
    seqan3::contrib::sdsl::load(this->sds, in);
    this->s_select = sd_vector<>::select_1_type(&sds);
    seqan3::contrib::sdsl::load(this->offset, in);
    seqan3::contrib::sdsl::load(this->span, in);
    in.close();
    return 0;
}

