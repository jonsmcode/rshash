#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/contrib/sdsl-lite.hpp>

#include "../source/minimiser_rev_hash_views.hpp"
#include "../source/minimiser_rev_hash_views2.hpp"


const uint8_t m_thres = 10;
const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
const uint64_t seed2 = 0x29'6D'BD'33'32'56'8C'64;


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
    const uint8_t k = 31;
    size_t N = 0;
    for(auto & record : input) {
        N += record.size();
    }
    // uint8_t m = std::log(N)/std::log(4);
    // if(m > 16)
    //     m = 16;
    // uint8_t m = 15;
    uint8_t m = 16;
    std::cout << "m: " << +m << '\n';

    // const size_t span = k-m;
    const size_t span = 1000000;

    auto view = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});

    const uint64_t M = 1ULL << (m+m-1);
    const uint64_t mask = compute_mask(m);

    std::cout << "scanning minimizers...\n";
    seqan3::contrib::sdsl::bit_vector r = seqan3::contrib::sdsl::bit_vector(M, 0);
    
    for(auto & record : input) {
        for(auto && minimiser : record | view) { // todo: simple minimiser view
            r[minimiser.minimiser_value] = 1;
        }
    }
    seqan3::contrib::sdsl::rank_support_v<1> r_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r);

    std::cout << "counting minimizers...\n";
    size_t c = r_rank(M);
    uint8_t* count = new uint8_t[c];
    std::memset(count, 0, c*sizeof(uint8_t));

    uint64_t kmers = 0;
    uint64_t n = 0;
    size_t longest_span = 0;
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            size_t i = r_rank(minimiser.minimiser_value);
            size_t o = minimiser.occurrences;

            size_t w = o/span + 1;
            if(count[i] + w >= m_thres) {
                count[i] = m_thres;
            }
            else {
                count[i] += w;
            }
            kmers += o;
            longest_span = std::max(o, longest_span);
            n += w;
        }
    }

    std::cout << "filling R_1 and R_2...\n";
    seqan3::contrib::sdsl::bit_vector r1 = seqan3::contrib::sdsl::bit_vector(M, 0);
    seqan3::contrib::sdsl::bit_vector r2 = seqan3::contrib::sdsl::bit_vector(M, 0);

    auto view12 = srindex::views::two_minimisers_and_window_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});

    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view12) {
            size_t i = r_rank(minimiser.minimiser1_value);
            if(count[i] < m_thres)
                r1[minimiser.minimiser1_value] = 1;
            else
                r2[minimiser.minimiser2_value] = 1;
        }
    }
    seqan3::contrib::sdsl::rank_support_v<1> r1_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r1);
    seqan3::contrib::sdsl::rank_support_v<1> r2_rank = seqan3::contrib::sdsl::rank_support_v<1>(&r2);

    seqan3::contrib::sdsl::util::clear(r);
    delete[] count;

    std::cout << "count minimizers1...\n";
    size_t c1 = r1_rank(M);
    uint8_t* count1 = new uint8_t[c1];
    std::memset(count1, 0, c1*sizeof(uint8_t));
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view) {
            if(r1[minimiser.minimiser_value]) {
                size_t i = r1_rank(minimiser.minimiser_value);
                count1[i] += minimiser.occurrences/span + 1;
            }
        }
    }
    uint64_t n1 = 0;
    for(uint64_t i=0; i < c1; i++)
        n1 += count1[i];

    std::cout << "count minimizers2...\n";
    size_t c2 = r2_rank(M);
    uint32_t* count2 = new uint32_t[c2];
    std::memset(count2, 0, c2*sizeof(uint32_t));

    auto view2 = srindex::views::minimiser_hash_and_positions({.minimiser_size = m, .window_size = k, .seed=seed2});
    for(auto & sequence : input) {
        for(auto && minimiser : sequence | view2) {
            if(r2[minimiser.minimiser_value]) {
                size_t i = r2_rank(minimiser.minimiser_value);
                count2[i] += minimiser.occurrences/span+1;
            }
        }
    }

    uint64_t n2 = 0;
    for(uint64_t i=0; i < c2; i++)
        n2 += count2[i];

    uint64_t thresholds[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    uint64_t* counter1 = new uint64_t[9];
    std::memset(counter1, 0, 9*sizeof(uint64_t));

    for(uint64_t i=0; i < c; i++) {
        if(count1[i] < m_thres) {
            for(int j=8; j >= 0; j--) {
                if(count1[i] >= thresholds[j]) {
                    counter1[j]++;
                    break;
                }
            }
        }
    }

    std::cout << "\nminimizer r1:\n";
    for(int j=0; j < 9; j++) {
        std::cout << "occurrences " << thresholds[j] << ": " << (double) counter1[j]/(c1+c2)*100 << "% abs: " << counter1[j] << "\n";
    }

    delete[] count1;
    delete[] counter1;

    uint64_t thresholds2[15] = {1, 2, 3, 4, 5, 6, 10, 20, 30, 40, 50, 100, 1000, 10000, 100000};
    uint64_t* counter2 = new uint64_t[15];
    std::memset(counter2, 0, 15*sizeof(uint64_t));
    for(uint64_t i=0; i < c2; i++) {
        for(int j=14; j >= 0; j--) {
            if(count2[i] >= thresholds2[j]) {
                counter2[j]++;
                break;
            }
        }
    }
    
    std::cout << "\nminimizer r2:\n";
    for(int j=0; j < 5; j++) {
        std::cout << "occurrences " << thresholds2[j] << ": " << (double) counter2[j]/(c1+c2)*100 << "% abs: " << counter2[j] << "\n";
    }
    for(int j=5; j < 14; j++) {
        std::cout << thresholds2[j] << "<= occurrences < " << thresholds2[j+1] << ": " << (double) counter2[j]/(c1+c2)*100 << "% abs: " << counter2[j] << "\n";
    }
    std::cout << "occurrences <= " << thresholds2[14] << ": " << (double) counter2[14]/(c1+c2)*100 << "% abs: " << counter2[14] << "\n";
    
    delete[] count2;
    delete[] counter2;

    std::cout << "\n====== report ======\n";
    std::cout << "text length: " << N << "\n";
    std::cout << "no kmers: " << kmers <<  '\n';
    std::cout << "no minimiser1: " << n1 << "  " << (double) n1/(n1+n2)*100 << "%\n";
    std::cout << "no distinct minimiser1: " << c1 << "  " << (double) c1/(c1+c2)*100 << "%\n";
    std::cout << "avg superkmers1: " << (double) n1/c1 <<  '\n';
    std::cout << "no minimiser2: " << n2 << "  " << (double) n2/(n1+n2)*100 << "%\n";
    std::cout << "no distinct minimiser2: " << c2 << "  " << (double) c2/(c1+c2)*100 << "%\n";
    std::cout << "avg superkmers2: " << (double) n2/c2 <<  '\n';
    std::cout << "density r1: " << (double) c1/M*100 << "%\n";
    std::cout << "density r2: " << (double) c2/M*100 << "%\n";
    std::cout << "longest span: " << longest_span <<  '\n';
    std::cout << "avg span: " << (double) kmers/n <<  '\n';

}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    stats(text);
}