#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

// #include "../source/minimiser_rev_hash_views.hpp"
#include "../source/minimiser_hash_views.hpp"


const uint8_t m_thres = 255;


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    uint64_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
        if(N >= 1000000000)
            break;
    }
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


void stats(const std::vector<std::vector<seqan3::dna4>> &input)
{
    uint8_t k = 31;
    size_t N = 0;
    for(auto & record : input) {
        N += record.size();
    }
    // uint8_t m = std::log(N)/std::log(4);
    // if(m > 16)
    //     m = 16;
    uint8_t m = 16;
    std::cout << "m: " << +m << '\n';

    auto view = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k});

    // const uint64_t M = 1ULL << (m+m-1);
    const uint64_t M = 1ULL << (m+m);
    const uint64_t mask = compute_mask(m);

    std::cout << "extracting minimizers...\n";
    seqan3::contrib::sdsl::bit_vector r1 = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r1[minimiser.minimiser_value] = 1;
        }
    }
    seqan3::contrib::sdsl::rank_support_v<1> r1_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r1);

    std::cout << "counting minimizers1...\n";
    size_t c = r1_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));
    std::unordered_map<uint64_t, uint32_t> cb;

    uint32_t max_occs = 0;
    uint64_t kmers = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r1_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            int w = o/k + 1;
            if(count[i] == m_thres) {
                cb[i] += w;
                max_occs = std::max(cb[i] + m_thres, max_occs);
            }
            else if(count[i] + w >= m_thres) {
                cb[i] = w - (m_thres - count[i]);
                count[i] = m_thres;
                r1[i] = 0;
            }
            else {
                count[i] += w;
            }
            kmers += o;
        }
    }

    uint64_t thresholds[10] = {1, 2, 3, 4, 5, 10, 100, 1000, 10000, 100000};
    uint64_t* counter = new uint64_t[10];
    std::memset(counter, 0, 10*sizeof(uint64_t));

    for(uint64_t i=0; i < c; i++) {
        if(count[i] == m_thres) {
            for(int j=9; j >= 0; j--) {
                if(cb[i]+m_thres >= thresholds[j]) {
                    counter[j]++;
                    break;
                }
            }
        }
        else {
            for(int j=9; j >= 0; j--) {
                if(count[i] >= thresholds[j]) {
                    counter[j]++;
                    break;
                }
            }
        }
        
    }

    delete[] count;

    std::cout << "====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    // std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "no distinct minimiser: " << c <<  '\n';
    // std::cout << "freq minimiser: " << cb.size() << ", " << (float) cb.size()/n*100 << "%\n";
    std::cout << "density r: " << (double) c/M*100 << "%\n";

    std::cout << "\n\nminimizer:\n";
    for(int j=0; j < 10; j++) {
        std::cout << "occurrences " << thresholds[j] << ": " << (double) counter[j]/c*100 << "% abs: " << counter[j] << "\n";
    }
    std::cout << "max. occurrences: " << max_occs << "\n";

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}