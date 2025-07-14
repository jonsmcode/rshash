#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../source/minimiser_rev_hash_views.hpp"


const uint8_t m_thres = 255;
// const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
// const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;
const uint64_t seed2 = 1;


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

static inline constexpr uint64_t compute_mask(uint64_t const size) {
        assert(size > 0u);
        assert(size <= 64u);

        if (size == 64u)
            return std::numeric_limits<uint64_t>::max();
        else
            return (uint64_t{1u} << (size)) - 1u;
    }


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


uint64_t load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    uint64_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
    }
    return N;
}


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    const uint8_t k = 31;
    size_t N = 0;
    for(auto & record : input) {
        N += record.size();
    }
    
    uint8_t m = 16;

    const size_t span = 100;

    auto view = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed2});

    const uint64_t M = 1ULL << (m+m);
    seqan3::contrib::sdsl::bit_vector r = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view)
            r[minimiser.minimiser_value] = 1;
    }
    seqan3::contrib::sdsl::rank_support_v<1> r_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r);
    
    size_t c = r_rank(M);
    uint32_t* count = new uint32_t[c];
    std::memset(count, 0, c*sizeof(uint32_t));

    uint32_t max_occs = 0;
    uint64_t max_minimizer;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            count[i] += minimiser.occurrences/span + 1;
            if(count[i] > max_occs) {
                max_occs = count[i];
                max_minimizer = minimiser.minimiser_value;
            }
        }
    }

    seqan3::debug_stream << "hashed minimizer: " << kmer_to_string(max_minimizer, m) << " appearing in " << max_occs << " super-kmers\n";

    std::unordered_map<uint64_t, uint32_t> freq_minimizers;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            if(!freq_minimizers.contains(minimiser.minimiser_value))
                freq_minimizers[minimiser.minimiser_value] = 0;
            freq_minimizers[minimiser.minimiser_value] += minimiser.occurrences/span + 1;
        }
    }

    std::vector<std::pair<uint64_t, uint32_t>> vec(freq_minimizers.begin(), freq_minimizers.end());
    std::sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) {return a.second > b.second;});

    auto view2 = srindex::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k, .seed=seed2});

    int n = 5;
    for (const auto& p : vec) {
        seqan3::debug_stream << "hashed minimizer: " << kmer_to_string(p.first, m) << " appearing in " << p.second << " super-kmers, covered by ";
        std::unordered_set<uint64_t> freq_kmers;
        for(auto & sequence : input) {
            for(auto && minimiser : sequence | view2) {
                if(minimiser.minimiser_value == p.first)
                    freq_kmers.insert(minimiser.window_value);
            }
        }
        seqan3::debug_stream << freq_kmers.size() << " kmers\n";
        if(n == 0)
            break;
        n--;
    }



    // for(auto & sequence : input) {
    //     for(auto && minimiser : sequence | view2) {
    //             size_t end = minimiser.range_position + k + minimiser.occurrences;
    //             if(end > sequence.size())
    //                 end = sequence.size();
    //             for(size_t j = minimiser.range_position; j < end; j++) {
    //                 seqan3::debug_stream << sequence[j];
    //             }
    //             uint64_t kmer;
    //             for(size_t j = minimiser.range_position; j < minimiser.range_position + k; j++) {
    //                 uint64_t const base = seqan3::to_rank(sequence[j]);
    //                 kmer = ((kmer << 2) | base);
    //             }
    //             freq_kmers.insert(kmer);
    //             for(size_t j = minimiser.range_position + k; j < end; j++) {
    //                 uint64_t const base = seqan3::to_rank(sequence[j]);
    //                 kmer = ((kmer << 2) | base) & kmer_mask;
    //                 freq_kmers.insert(kmer);
    //             }
    //             seqan3::debug_stream << '\n';

    //     }
    // }

    // for (const auto& p : vec) {
    //     seqan3::debug_stream << "hashed minimizer: " << kmer_to_string(p.first, m) << " appearing in " << p.second << " super-kmers, covered by " << freq_kmers[p.first].size() << " kmers\n";
    // }

    // for(auto & sequence : input) {
    //     for(auto && minimiser : sequence | view) {
    //         size_t i = r_rank(minimiser.minimiser_value);
    //         if(count[i] == m_thres) {
    //             size_t end = minimiser.range_position + k + minimiser.occurrences;
    //             if(end > sequence.size())
    //                 end = sequence.size();
    //             for(uint64_t j = minimiser.range_position; j < end; j++) {
    //                 seqan3::debug_stream << sequence[j];
    //             }
    //             seqan3::debug_stream << '\n';
    //             seqan3::debug_stream << kmer_to_string(minimiser.minimiser_value, m) << ' ' << m_thres + cb[i] << '\n';
    //         }
    //     }
    // }
    // seqan3::debug_stream << '\n';

    // delete[] count;

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/bigdata/ag_abi/jonas/kmerdict/skmers2.fasta";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}