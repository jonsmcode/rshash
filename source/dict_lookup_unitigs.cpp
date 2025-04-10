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
    // todo: simons bitvector
    r = bit_vector(M);
    // todo: faster minimiser only view

    size_t longest_sequence = 0;
    size_t N = 0;
    uint64_t n_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r[minimiser.minimiser_value] = 1;
        }
        longest_sequence = std::max(longest_sequence, record.size());
        n_sequences++;
        N += record.size();
    }
    r_rank = rank_support_v<1>(&r);

    size_t c = r_rank(M);

    std::cout << "counting minimizers...\n";
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

    size_t sequences_width = std::bit_width(n_sequences);
    size_t positions_width = std::bit_width(longest_sequence);
    sequences.width(sequences_width);
    sequences.resize(n);
    positions.width(positions_width);
    positions.resize(n);

    std::cout << "allocated " << n*(sequences_width + positions_width)/8 << " bytes for index\n";

    std::memset(count, 0, c*sizeof(uint8_t));
    // faster?
    for (const auto& [key, value] : cb)
        cb[key] = 0;

    std::cout << "building index...\n";

    uint64_t sequence_id = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view) {
            uint64_t i = r_rank(minimiser.minimiser_value);
            uint64_t s = simple_select.select(i);
            int o = minimiser.occurrences;
            int j = 0;
            while(o > k) {
                if(count[i] == 255) {
                    positions[s + count[i] + cb[i]] = minimiser.range_position + j*k;
                    sequences[s + count[i] + cb[i]] = sequence_id;
                    cb[i]++;
                }
                else {
                    positions[s + count[i]] = minimiser.range_position + j*k;
                    sequences[s + count[i]] = sequence_id;
                    count[i]++;
                }
                o -= k;
                j++;
            }
            if(count[i] == 255) {
                positions[s + count[i] + cb[i]] = minimiser.range_position + j*k;
                sequences[s + count[i] + cb[i]] = sequence_id;
                cb[i]++;
            }
            else {
                positions[s + count[i]] = minimiser.range_position + j*k;
                sequences[s + count[i]] = sequence_id;
                count[i]++;
            }
        }
        sequence_id++;
    }

    delete[] count; // delete cb?


    // for(auto & record : input) {
    //     std::ranges::move(record, std::back_inserter(text));
    // }
    text = input;

    // report
    std::cout << "====== report ======\n";
    std::cout << "number sequences: " << n_sequences << "\n";
    std::cout << "longest sequence: " << longest_sequence << "\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "avg sequence length: " << (float) N/n_sequences << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no distinct kmers: " << c << '\n';
    std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "freq minimiser: " << cb.size() << '\n';
    std::cout << "density r: " << (float) c/M*100 << "%\n";
    std::cout << "density s: " << (float) simple_select.bitCount()/(n+1)*100 <<  "%\n";
    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "sequences: " << (float) n*sequences_width/kmers << "\n";
    std::cout << "positions: " << (float) n*positions_width/kmers << "\n";
    std::cout << "text: " << (float) 2*N/kmers << "\n";
    std::cout << "r: " << (float) M/kmers << "\n";
    std::cout << "s: " << (float) (n+1)/kmers << "\n";
    std::cout << "total: " << (float) (n*(sequences_width+positions_width)+2*N+M+n+1)/kmers << "\n";

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
                size_t seq = sequences[p+i];
                size_t pos = positions[p+i];
                for (int j=pos; j < pos+k; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[seq][j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                }
                if(minimiser.window_value == hash) {
                    occurences++;
                    break;
                }
                for(int j=pos+k; j < std::min(pos+2*k, text[seq].size()); j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[seq][j]);
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


int LookupUnitigsDictionary::streaming_query(const std::vector<seqan3::dna4> &query,
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> &result)
{
    auto query_view = bsc::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k});

    uint64_t kmer = 0;
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
                size_t seq = sequences[p+i];
                size_t pos = positions[p+i];
                for (int j=pos; j < pos+k; j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[seq][j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                }
                if(minimiser.window_value == hash) {
                    result.push_back({kmer, seq, pos});
                    break;
                }
                for(int j=pos+k; j < std::min(pos+2*k, text[seq].size()); j++) {
                    uint64_t const new_rank = seqan3::to_rank(text[seq][j]);
                    hash <<= 2;
                    hash |= new_rank;
                    hash &= mask;
                    if(minimiser.window_value == hash) {
                        found = true;
                        result.push_back({kmer, seq, j-k+1});
                        break;
                    }
                }
            }

        }
        // todo: else ht lookup

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
    seqan3::contrib::sdsl::serialize(this->sequences, out);
    seqan3::contrib::sdsl::serialize(this->positions, out);

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
    seqan3::contrib::sdsl::load(this->sequences, in);
    seqan3::contrib::sdsl::load(this->positions, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);

    in.close();
    return 0;
}

