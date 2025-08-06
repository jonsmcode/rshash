#include <filesystem>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <cereal/archives/binary.hpp>

#include "rsindex.hpp"
#include "io.hpp"
#include "minimiser_rev_xor_views.hpp"


const uint64_t seed = 0x8F'3F'73'B5'CF'1C'9A'DE;


static inline constexpr uint64_t compute_mask(uint64_t const size)
{
    assert(size > 0u);
    assert(size <= 64u);

    if (size == 64u)
        return std::numeric_limits<uint64_t>::max();
    else
        return (uint64_t{1u} << size) - 1u;
}


RSIndexComp::RSIndexComp() {}

RSIndexComp::RSIndexComp(uint8_t const k, uint8_t const m, uint8_t const m_thres) {
    this->k = k;
    this->m = m;
    this->m_thres = m_thres;
}


int RSIndexComp::build(const std::vector<std::vector<seqan3::dna4>> &input)
{
    // auto view1 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed});
    auto view1 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed});

    // assert(m <= 32);
    const uint64_t M = 1ULL << (m+m);

    std::cout << "find minimizers...\n";
    bit_vector rtmp = bit_vector(M, 0);

    size_t N = 0;
    uint64_t no_sequences = 0;
    for(auto & record : input) {
        for(auto && minimiser : record | view1) {
            rtmp[minimiser.minimiser_value] = 1;
        }
        N += record.size();
        no_sequences++;
    }
    rank_support_v<1> rtmp_rank = rank_support_v<1>(&rtmp);

    std::cout << "get sequences...\n";
    bit_vector sequences = bit_vector(N+1, 0);
    size_t j = 0;
    sequences[0] = 1;
    for(uint64_t i=0; i < no_sequences; i++) {
        j += input[i].size();
        sequences[j] = 1;
    }
    endpoints = seqan3::contrib::sdsl::sd_vector<>(sequences);
    endpoints_rank = seqan3::contrib::sdsl::rank_support_sd<>(&endpoints);
    endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);

    std::cout << "count minimizers...\n";
    size_t ctmp = rtmp_rank(M);
    uint8_t* counttmp = new uint8_t[ctmp];
    std::memset(counttmp, 0, ctmp*sizeof(uint8_t));

    uint64_t kmers = 0;
    uint64_t skmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            size_t i = rtmp_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;
            size_t w = o/span + 1;
            counttmp[i] += w;
            if(counttmp[i] > m_thres)
                counttmp[i] = m_thres;
            kmers += o;
            skmers += w;
        }
    }

    std::cout << "filling R...\n";
    bit_vector r_ = bit_vector(M, 0);
    for(auto & sequence : input) {
        for(auto && minimisers : sequence | view1) {
            if(counttmp[rtmp_rank(minimisers.minimiser_value)] < m_thres)
                r_[minimisers.minimiser_value] = 1;
        }
    }
    r = sd_vector<>(r_);
    r_rank = rank_support_sd<>(&r);

    delete[] counttmp;

    std::cout << "count minimisers again...\n";
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view1) {
            if(r[minimiser.minimiser_value]) {
                size_t i = r_rank(minimiser.minimiser_value);
                count[i] += minimiser.occurrences/span + 1;
            }
        }
    }
    uint64_t n = 0;
    for(size_t i=0; i < c; i++)
        n += count[i];

    std::cout << "filling bitvector S...\n";
    s = bit_vector(n+1, 0);
    s[0] = 1;
    j = 0;
    for (size_t i=0; i < c; i++) {
        j += count[i];
        s[j] = 1;
    }
    s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), n+1, 3);

    std::cout << "filling offsets...\n";
    const size_t offset_width = std::bit_width(N);
    offsets.width(offset_width);
    offsets.resize(n);

    std::memset(count, 0, c*sizeof(uint8_t));

    size_t length = 0;
    for(auto & sequence : input) {
        for (auto && minimiser : sequence | view1) {
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

    delete[] count;

    std::cout << "get frequent skmers...\n";
    std::vector<std::vector<seqan3::dna4>> freq_skmers;
    for(auto & sequence : input) {
        size_t start_position;
        bool level_up;
        bool current_level_up;

        for(auto && minimiser : sequence | view1) {
            level_up = r[minimiser.minimiser_value];
            break;
        }
        if(!level_up)
            start_position = 0;
        for(auto && minimiser : sequence | view1) {
            current_level_up = r[minimiser.minimiser_value];

            if(level_up && !current_level_up)
                start_position = minimiser.range_position;
            if(!level_up && current_level_up) {
                std::vector<seqan3::dna4> skmer;
                for(size_t i=start_position; i < minimiser.range_position+k; i++)
                    skmer.push_back(sequence[i]);
                freq_skmers.push_back(skmer);
            }

            level_up = current_level_up;
        }
        if(!current_level_up) {
            std::vector<seqan3::dna4> skmer;
            for(size_t i=start_position; i < sequence.size(); i++)
                skmer.push_back(sequence[i]);
            freq_skmers.push_back(skmer);
        }
    }
    
    size_t len_rem_seqs = 0;
    for(auto & sequence : freq_skmers)
        len_rem_seqs += sequence.size();
    std::cout << "remaining superkmers " << freq_skmers.size() << " (" << (double) freq_skmers.size()/n*100 << "%) ";
    std::cout << "total length: " << len_rem_seqs << " (" << (double) len_rem_seqs/N*100 << "%)\n";

    std::cout << "build level 2, HT...\n";

    std::cout << "filling HT...\n";
    auto view3 = srindex::views::xor_minimiser_and_window({.minimiser_size = m, .window_size = k, .seed=seed});
    for(auto & sequence : freq_skmers) {
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
    std::cout << "no skmers: " << skmers <<  '\n';
    
    std::cout << "no minimiser: " << n << "\n";
    std::cout << "no distinct minimiser: " << ctmp << "\n";
    std::cout << "minimiser going to HT: " << ctmp-c << "  " << (double) (ctmp-c)/ctmp*100 << "%\n";
    std::cout << "avg superkmers: " << (double) n/c <<  '\n';
    // std::cout << "no minimiser HT: " << hashmap.size() << " " << (double)hashmap.size()/n*100 << "%\n";
    std::cout << "no kmers HT: " << hashmap.size() << " " << (double) hashmap.size()/kmers*100 << "%\n";
    std::cout << "density r1: " << (double) c/M*100 << "%\n";
    std::cout << "density s1: " << (double) s_select.bitCount()/(n+1)*100 <<  "%\n";

    std::cout << "\nspace per kmer in bit:\n";
    std::cout << "text: " << (double) 2*N/kmers << "\n";
    std::cout << "endpoints: " << (double) 8*size_in_bytes(endpoints)/kmers << "\n";
    std::cout << "offsets: " << (double) n*offset_width/kmers << "\n";
    std::cout << "Hashtable: " << (double) 64*hashmap.size()/kmers << "\n";
    std::cout << "R: " << (double) 8*size_in_bytes(r)/kmers << "\n";
    std::cout << "S: " << (double) (n+1)/kmers << "\n";

    std::cout << "total: " << (double) (n*offset_width+2*N+8*size_in_bytes(r)+n+1+8*size_in_bytes(endpoints)+64*hashmap.size())/kmers << "\n";

    return 0;
}


inline void RSIndexComp::fill_buffer(std::vector<uint64_t> &buffer, const uint64_t mask, size_t p, size_t q)
{
    for(uint64_t i = 0; i < q-p; i++) {
        uint64_t hash = 0;
        size_t o = offsets[p+i];
        size_t next_endpoint = endpoints_select(endpoints_rank(o+1)+1);
        size_t e = o+k+span;
        if(e > next_endpoint)
            e = next_endpoint;
        for (uint64_t j=o; j < o+k; j++) {
            uint64_t const new_rank = seqan3::to_rank(text[j]);
            hash <<= 2;
            hash |= new_rank;
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


inline bool lookup_serial(std::vector<uint64_t> &array, uint64_t query, uint64_t queryrc, size_t &last_found) {
    for(size_t i=last_found+1; i < array.size(); i++) {
        if(array[i] == query || array[i] == queryrc) {
            last_found = i;
            return true;
        }
    }
    for(size_t i=0; i < last_found+1; i++) {
        if(array[i] == query || array[i] == queryrc) {
            last_found = i;
            return true;
        }
    }
    return false;
}


uint64_t RSIndexComp::streaming_query(const std::vector<seqan3::dna4> &query)
{
    auto view = srindex::views::xor_minimiser_and_window({.minimiser_size = m, .window_size = k, .seed=seed});
    uint64_t occurences = 0;
    const uint64_t mask = compute_mask(2u * k);
    uint64_t current_minimiser=std::numeric_limits<uint64_t>::max();
    std::vector<uint64_t> buffer;
    size_t last_found = 0;

    for(auto && minimisers : query | view)
    {
        if(minimisers.minimiser_value == current_minimiser)
            occurences += lookup_serial(buffer, minimisers.window_value, minimisers.window_value_rev, last_found);
        else if(r[minimisers.minimiser_value]) {
            size_t minimizer_id = r_rank(minimisers.minimiser_value);
            size_t p = s_select.select(minimizer_id);
            size_t q = s_select.select(minimizer_id+1);

            buffer.clear();
            fill_buffer(buffer, mask, p, q);
            last_found = 0;
            occurences += lookup_serial(buffer, minimisers.window_value, minimisers.window_value_rev, last_found);
            current_minimiser = minimisers.minimiser_value;
        }
        else
            occurences += hashmap.contains(std::min<uint64_t>(minimisers.window_value, minimisers.window_value_rev));
    }

    return occurences;
}


int RSIndexComp::save(const std::filesystem::path &filepath) {
    std::ofstream out(filepath, std::ios::binary);
    seqan3::contrib::sdsl::serialize(this->k, out);
    seqan3::contrib::sdsl::serialize(this->m, out);
    seqan3::contrib::sdsl::serialize(r, out);
    seqan3::contrib::sdsl::serialize(s, out);
    seqan3::contrib::sdsl::serialize(this->offsets, out);
    seqan3::contrib::sdsl::serialize(this->endpoints, out);

    cereal::BinaryOutputArchive archive(out);
    archive(this->text);
    archive(this->hashmap);

    out.close();
    return 0;
}

int RSIndexComp::load(const std::filesystem::path &filepath) {
    std::ifstream in(filepath, std::ios::binary);
    seqan3::contrib::sdsl::load(this->k, in);
    seqan3::contrib::sdsl::load(this->m, in);
    seqan3::contrib::sdsl::load(r, in);
    seqan3::contrib::sdsl::load(s, in);
    seqan3::contrib::sdsl::load(this->offsets, in);
    seqan3::contrib::sdsl::load(this->endpoints, in);

    cereal::BinaryInputArchive archive(in);
    archive(this->text);
    archive(this->hashmap);

    in.close();

    r_rank = rank_support_sd<>(&r);
    this->s_select = sux::bits::SimpleSelect(reinterpret_cast<uint64_t*>(s.data()), s.size(), 3);
    endpoints_rank = rank_support_sd<>(&endpoints);
    endpoints_select = seqan3::contrib::sdsl::select_support_sd<>(&endpoints);
    
    return 0;
}

