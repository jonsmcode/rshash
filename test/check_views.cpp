#include <filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "../source/minimiser_rev_xor_hash_views.hpp"


const uint64_t seed1 = 0x8F'3F'73'B5'CF'1C'9A'DE;
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
    uint64_t n = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
        if(N >= 10000000)
            break;
        n++;
    }
    std::cout << n << '\n';
    return N;
}


void check(const std::vector<std::vector<seqan3::dna4>> &input)
{
    const uint8_t k = 31;
    uint8_t m = 16;
    std::cout << "m: " << +m << '\n';

    const size_t span = 100;

    auto view1 = srindex::views::two_minimisers_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});
    auto view2 = srindex::views::two_minimisers_and_window_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});
    auto view3 = srindex::views::two_minimisers_and_occurence_hash({.minimiser_size = m, .window_size = k, .seed1=seed1, .seed2=seed2});
    auto view4 = srindex::views::xor_minimiser_and_positions({.minimiser_size = m, .window_size = k, .seed=seed1});

    std::vector<uint64_t> minimisers1_1;
    std::vector<uint64_t> minimisers1_2;
    std::vector<uint64_t> minimisers2_1;
    std::vector<uint64_t> minimisers2_2;
    std::vector<uint64_t> minimisers3_1;
    std::vector<uint64_t> minimisers3_2;
    // std::vector<uint64_t> minimisers4;


    std::cout << "scanning minimizers 1...\n";

    for(auto & record : input) {
        for(auto && minimiser : record | view1) {
            minimisers1_1.push_back(minimiser.minimiser1_value);
            minimisers1_2.push_back(minimiser.minimiser2_value);
        }
    }

    std::cout << "scanning minimizers 2...\n";
    for(auto & record : input) {
        for(auto && minimiser : record | view2) {
            minimisers2_1.push_back(minimiser.minimiser1_value);
            minimisers2_2.push_back(minimiser.minimiser2_value);
        }
    }

    std::cout << "scanning minimizers 3...\n";

    for(auto & sequence : input) {
        uint64_t current_m1;
        uint64_t current_m2;
        for(auto && minimiser : sequence | view2) {
            current_m1 = minimiser.minimiser1_value;
            current_m2 = minimiser.minimiser2_value;
            break;
        }
        for(auto && minimiser : sequence | view2) {
            minimisers3_1.push_back(current_m1);
            minimisers3_2.push_back(current_m2);
            current_m1 = minimiser.minimiser1_value;
            current_m2 = minimiser.minimiser2_value;
        }
    }

    // std::cout << "scanning minimizers 3...\n";
    // for(auto & record : input) {
    //     for(auto && minimiser : record | view3) {
    //         minimisers3_1.push_back(minimiser.minimiser1_value);
    //         minimisers3_2.push_back(minimiser.minimiser2_value);
    //     }
    // }
    // std::cout << "scanning minimizers 4...\n";
    // for(auto & record : input) {
    //     for(auto && minimiser : record | view4) {
    //         minimisers4.push_back(minimiser.minimiser_value);
    //     }
    // }

    std::cout << "checking view 1 and 2...\n";

    if(minimisers1_1.size() != minimisers2_1.size()) {
        std::cout << "minimisers 1 size differ\n";
        return;
    }

    if(minimisers1_2.size() != minimisers2_2.size()) {
        std::cout << "minimisers 2 size differ\n";
        return;
    }

    for(size_t i=0; i < minimisers1_1.size(); i++)
        if(minimisers1_1[i] != minimisers2_1[i]) {
            std::cout << "minimisers 1 differ\n";
            return;
        }
    for(size_t i=0; i < minimisers1_2.size(); i++) {
        if(minimisers1_2[i] != minimisers2_2[i]) {
            std::cout << "minimisers 2 differ\n";
            return;
        }
    }

    std::cout << "checking view 1 and 3...\n";

    // for(size_t i=0; i < minimisers3_1.size(); i++) {
    //         // seqan3::debug_stream << kmer_to_string(minimisers1_1[i], m) << ' '<< kmer_to_string(minimisers3_1[i], m) << '\n';
    //     std::cout << minimisers1_1[i] << ' '<< minimisers3_1[i] << '\n';
    // }
    // for(size_t i=0; i < minimisers3_1.size(); i++) {
    //     if(minimisers1_1[i] != minimisers3_1[i]) {
    //         std::cout << i << ' ';
    //         // seqan3::debug_stream << i << ' ' << kmer_to_string(minimisers1_1[i], m) << ' '<< kmer_to_string(minimisers3_1[i], m) << '\n';
    //     }
    // }

    if(minimisers1_1.size() != minimisers3_1.size()) {
        std::cout << "minimisers 1 size differ " << minimisers1_1.size() << ' ' << minimisers3_1.size() << '\n';
        // return;
    }
    if(minimisers1_2.size() != minimisers3_2.size()) {
        std::cout << "minimisers 2 size differ " << minimisers1_2.size() << ' ' << minimisers3_2.size() << '\n';
        // return;
    }

    // for(size_t i=0; i < 100; i++)
    //     std::cout << i << " " << minimisers1_1[i] << ' ' << minimisers1_2[i] << "   " << minimisers3_1[i] << ' ' << minimisers3_2[i] << '\n';
    // for(size_t i=minimisers1_1.size()-10; i < minimisers1_1.size(); i++)
    //     std::cout << i << " " << minimisers1_1[i] << '\n';
    // for(size_t i=minimisers3_1.size()-10; i < minimisers3_1.size(); i++)
    //     std::cout << i << " " << minimisers3_1[i] << '\n';

    // for(size_t i=0; i < minimisers1_1.size(); i++)
    //     if(minimisers1_1[i] != minimisers3_1[i]) {
    //         std::cout << i << " ";
    //         // return;
    //         break;
    //     }
    // for(size_t i=0; i < minimisers1_2.size(); i++) {
    //     if(minimisers1_2[i] != minimisers3_2[i]) {
    //         std::cout << i << " ";
    //         // return;
    //         break;
    //     }
    // }

    // std::cout << "checking view 1 and 4...\n";

    // if(minimisers1_1.size() != minimisers4.size()) {
    //     std::cout << "minimisers 1 size differ " << minimisers1_1.size() << ' ' << minimisers4.size() << '\n';
    //     return;
    // }

    // for(size_t i=0; i < minimisers1_1.size(); i++)
    //     if(minimisers1_1[i] != minimisers4[i]) {
    //         std::cout << "minimisers 1 differ\n";
    //         return;
    //     }
}



int main(int argc, char** argv)
{
    std::filesystem::path path = "/Users/adm_js4718fu/datasets/unitigs/human.k31.unitigs.fa.ust.fa.gz";
    std::vector<std::vector<seqan3::dna4>> text;
    load_file(path, text);
    check(text);
}