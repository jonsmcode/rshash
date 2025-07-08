#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

#include "../source/minimiser_rev_hash_views.hpp"
// #include "../source/minimiser_hash_views.hpp"


const uint8_t m_thres = 10;
const uint64_t seed = 0x8F'3F'73'B5'CF'1C'9A'DE;


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

    auto view = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed});

    const uint64_t M = 1ULL << (m+m);

    std::cout << "extracting minimizers...\n";
    seqan3::contrib::sdsl::bit_vector r = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view) {
            r[minimiser.minimiser_value] = 1;
        }
    }
    seqan3::contrib::sdsl::rank_support_v<1> r_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r);

    std::cout << "counting minimizers...\n";
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));
    std::unordered_map<uint64_t, uint32_t> cb;

    uint32_t max_occs = 0;
    uint64_t kmers = 0;
    uint64_t n = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            const uint64_t w = o/k + 1;
            if(count[i] == m_thres) {
                cb[i] += w;
                max_occs = std::max(cb[i] + m_thres, max_occs);
            }
            else if(count[i] + w >= m_thres) {
                cb[i] = w - (m_thres - count[i]);
                count[i] = m_thres;
            }
            else {
                count[i] += w;
            }
            kmers += o;
            n += w;
        }
    }

    std::cout << "counting freq kmers...\n";

    std::unordered_set<uint64_t> freq_kmers;
    auto view2 = srindex::views::minimiser_and_window_hash({.minimiser_size = m, .window_size = k, .seed=seed});
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view2) {
            size_t i = r_rank(minimiser.minimiser_value);
            if(count[i] == m_thres) {
                freq_kmers.insert(minimiser.window_value);
            }
        }
    }


    uint64_t thresholds[11] = {1, 2, 3, 4, 5, 6, 10, 100, 1000, 10000, 100000};
    uint64_t* counter = new uint64_t[11];
    std::memset(counter, 0, 11*sizeof(uint64_t));

    for(uint64_t i=0; i < c; i++) {
        if(count[i] == m_thres) {
            for(int j=10; j >= 0; j--) {
                if(cb[i]+m_thres >= thresholds[j]) {
                    counter[j]++;
                    break;
                }
            }
        }
        else {
            for(int j=10; j >= 0; j--) {
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
    std::cout << "no minimiser: " << n <<  '\n';
    std::cout << "no distinct minimiser: " << c <<  '\n';
    std::cout << "avg superkmers: " << (double) n/c <<  '\n';
    std::cout << "density r: " << (double) c/M*100 << "%\n";

    std::cout << "\n\nminimizer:\n";
    for(int j=0; j < 5; j++) {
        std::cout << "occurrences " << thresholds[j] << ": " << (double) counter[j]/c*100 << "% abs: " << counter[j] << "\n";
    }
    for(int j=5; j < 10; j++) {
        std::cout << thresholds[j] << "<= occurrences < " << thresholds[j+1] << ": " << (double) counter[j]/c*100 << "% abs: " << counter[j] << "\n";
    }
    std::cout << "occurrences <= " << thresholds[10] << ": " << (double) counter[10]/c*100 << "% abs: " << counter[10] << "\n";
    std::cout << "max. occurrences: " << max_occs << "\n";

    std::cout << "\n\nfreq kmers: " << freq_kmers.size() << " " << (double) freq_kmers.size()/kmers*100 << "%\n";

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}