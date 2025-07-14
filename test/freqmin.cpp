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
    std::cout << "m: " << +m << '\n';

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
    
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));
    std::unordered_map<uint64_t, uint32_t> cb;
    uint64_t max_occs = 0;
    uint64_t max_minimiser = 0;

    uint64_t kmers = 0;
    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            const uint64_t w = o/span + 1;
            if(count[i] == m_thres) {
                cb[i] += w;
                max_occs = std::max<uint64_t>(cb[i], max_occs);
                max_minimiser = minimiser.minimiser_value;
            }
            else if(count[i] + w >= m_thres) {
                cb[i] = w - (m_thres - count[i]);
                count[i] = m_thres;
            }
            else
                count[i] += w;
            kmers += o;
            n += w;
        }
    }

    seqan3::debug_stream << kmer_to_string(max_minimiser, m) << " occs: " << max_occs << "\n\n";
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            if(minimiser.minimiser_value == max_minimiser) {
                size_t end = minimiser.range_position + k + minimiser.occurrences;
                if(end > sequence.size())
                    end = sequence.size();
                for(size_t j = minimiser.range_position; j < end; j++) {
                    seqan3::debug_stream << sequence[j];
                }
                seqan3::debug_stream << '\n';
            }
        }
    }
    seqan3::debug_stream << '\n';

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

    delete[] count;

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/bigdata/ag_abi/jonas/kmerdict/skmers2.fasta";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}