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

    std::cout << "extracting minimizers...\n";
    // todo: test simons bitvector
    r = bit_vector(M, 0);

    size_t N = 0;
    uint64_t number_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r[minimiser.minimiser_value] = 1;
        }
        N += record.size();
        number_sequences++;
    }
    r_rank = rank_support_v<1>(&r);

    std::cout << "filling bitvector sequences and sequence lengths...\n";
    // todo: Elias Fano
    sequences = bit_vector(N, 0);
    endpoints.width(std::bit_width(N));
    endpoints.resize(number_sequences+1);
    endpoints[0] = 0;
    uint64_t j = 0;
    for(int i=0; i < number_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
        endpoints[i+1] = j;
    }
    sequences_rank = rank_support_v<1>(&sequences);

    std::cout << "counting minimizers...\n";
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    // todo: init hashtable with adequate size, better ht?
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
    std::cout << "filling bitvector S...\n";
    s = bit_vector(n+1, 0);
    s[0] = 1;
    j = 0;
    for (size_t i=0; i < c; i++) {
        j += count[i];
        if(count[i] == 255)
            j += cb[i];
        s[j] = 1;
    }
    s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), n+1, 3);

    size_t offset_width = std::bit_width(N);
    offsets.width(offset_width);
    offsets.resize(n);

    std::cout << "allocated " << n*(offset_width)/8000000 << " Mega bytes for O\n";

    std::memset(count, 0, c*sizeof(uint8_t));
    // faster?
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    std::cout << "filling O...\n";

    uint64_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view) {
            uint64_t i = r_rank(minimiser.minimiser_value);
            uint64_t s = s_select.select(i);
            int o = minimiser.occurrences;
            int j = 0;
            while(o > k) {
                if(count[i] == 255) {
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
            if(count[i] == 255) {
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

    delete[] count; // delete cb?


    for(auto & record : input) {
        std::ranges::move(record, std::back_inserter(text));
    }

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no distinct kmers: " << c << '\n';
    std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "freq minimiser: " << cb.size() << '\n';
    std::cout << "density r: " << (float) c/M*100 << "%\n";
    std::cout << "density s: " << (float) s_select.bitCount()/(n+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "offsets: " << (float) n*offset_width/kmers << "\n";
    std::cout << "text: " << (float) 2*N/kmers << "\n";
    std::cout << "R: " << (float) M/kmers << "\n";
    std::cout << "S: " << (float) (n+1)/kmers << "\n";
    std::cout << "Sequences: " << (float) N/kmers << "\n";
    std::cout << "endpoints: " << (float) (std::bit_width(N)*number_sequences)/kmers << "\n";
    std::cout << "total: " << (float) (n*offset_width+2*N+M+n+1+N+std::bit_width(N)*number_sequences)/kmers << "\n";

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
            size_t minimizer_id = r_rank(minimiser.minimiser_value);
            size_t p = s_select.select(minimizer_id);
            size_t q = s_select.select(minimizer_id+1);
            size_t b = q - p;

            bool found = false;
            for(uint64_t i = 0; !found && i < b; i++) {
                uint64_t hash = 0;
                size_t o = offsets[p+i];
                size_t sequence_id = sequences_rank(o);
                for (uint64_t j=o; j < o+k; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                }
                if(minimiser.window_value == hash) {
                    occurences++;
                    break;
                }
                for(uint64_t j=o+k; j < std::min<uint64_t>(o + 2*k, (size_t) endpoints[sequence_id+1]); j++) {
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
    }

    return occurences;
}


int LookupUnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &query,
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &result)
{
    auto view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

    const uint64_t mask = compute_mask(k);
    uint64_t kmer = 0;

    for(auto && minimiser : query | view) {
        if(r[minimiser.minimiser_value]) {
            uint64_t minimizer_id = r_rank(minimiser.minimiser_value);
            uint64_t p = s_select.select(minimizer_id);
            uint64_t q = s_select.select(minimizer_id+1);

            bool found = false;
            for(uint64_t i = 0; !found && i < q-p; i++) {
                uint64_t hash = 0;
                uint64_t o = offsets[p+i];
                uint64_t sequence_id = sequences_rank(o);
                for (uint64_t j=o; j < o+k; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                }
                if(minimiser.window_value == hash) {
                    result.push_back({kmer, sequence_id, o-endpoints[sequence_id]});
                    break;
                }
                uint64_t e = o+k+k;
                if(e > endpoints[sequence_id+1])
                    e = endpoints[sequence_id+1];
                for(uint64_t j=o+k; j < e; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                    if(minimiser.window_value == hash) {
                        found = true;
                        result.push_back({kmer, sequence_id, j-k+1-endpoints[sequence_id]});
                        break;
                    }
                }
            }

        }
        // else {
        //     cb.find(minimiser.window_value)
        // }
        kmer++;
    }

    return 0;
}


int LookupUnitigsDictionary::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(r, out); // save r_rank?
    seqan3::contrib::sdsl::serialize(s, out);
    seqan3::contrib::sdsl::serialize(sequences, out);
    seqan3::contrib::sdsl::serialize(this->offsets, out);
    seqan3::contrib::sdsl::serialize(this->endpoints, out);

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
    this->s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), s.size(), 3);
    seqan3::contrib::sdsl::load(sequences, in);
    sequences_rank = rank_support_v<1>(&sequences);
    seqan3::contrib::sdsl::load(this->offsets, in);
    seqan3::contrib::sdsl::load(this->endpoints, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);

    in.close();
    return 0;
}

