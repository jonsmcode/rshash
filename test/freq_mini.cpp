#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <sharg/all.hpp>

#include "../source/minimiser_views.hpp"

using namespace seqan3::contrib::sdsl;

struct cmd_arguments {
    std::filesystem::path i{};
    uint8_t k{31};
    uint8_t m{17};
};

void initialise_argument_parser(sharg::parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, sharg::config{.short_id = 'i', .long_id = "input", .description = "provide input file"});
    parser.add_option(args.k, sharg::config{.short_id = 'k', .long_id = "kmer", .description = "k-mer length"});
    parser.add_option(args.m, sharg::config{.short_id = 'm', .long_id = "mini", .description = "minimiser length"});
}


struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath,
               std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &output)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        seqan3::bitpacked_sequence<seqan3::dna4> seq;
        seq.assign(record.sequence().begin(), record.sequence().end());
        output.push_back(std::move(seq));
    }
}


std::string decode_kmer(uint64_t kmer, size_t k)
{
    std::string result(k, 'A');

    for (size_t i = 0; i < k; ++i)
    {
        uint8_t bits = kmer & 0b11;

        switch (bits)
        {
            case 0b00: result[k - 1 - i] = 'A'; break;
            case 0b01: result[k - 1 - i] = 'C'; break;
            case 0b10: result[k - 1 - i] = 'G'; break;
            case 0b11: result[k - 1 - i] = 'T'; break;
        }

        kmer >>= 2;
    }

    return result;
}


uint64_t count_minimizers(const std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> &sequences,
    const uint8_t minimiser_length, const uint8_t k, const uint64_t seed, uint64_t &max_minimizer, uint64_t &kmer)
{
    // std::cout << "computing minimizers...\n";
    std::vector<uint64_t> minimizers;
    auto minimiserview = rshash::views::xor_minimiser_and_positions({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    for(auto & sequence : sequences) {
        for(auto && minimiser : sequence | minimiserview)
            minimizers.emplace_back(minimiser.minimiser_value);
    }

    // std::cout << "sorting minimizers...\n";
    std::sort(minimizers.begin(), minimizers.end());

    // std::cout << "counting minimizers...\n";
    uint64_t current_minimizer = minimizers[0];
    uint64_t max_occurences = 0;
    uint64_t occurences = 0;
    for(uint64_t minimizer : minimizers) {
        if(minimizer != current_minimizer) {
            if(occurences > max_occurences) {
                max_occurences = occurences;
                max_minimizer = current_minimizer;
            }
            current_minimizer = minimizer;
            occurences = 1;
        }
        else
            occurences++;
    }
    if(minimizers.back() == current_minimizer && occurences > max_occurences) {
        max_occurences = occurences;
        max_minimizer = current_minimizer;
    }

    minimizers.clear();

    // auto minimiser_kmer_view = rshash::views::xor_minimiser_and_window({.minimiser_size = minimiser_length, .window_size = k, .seed=seed});
    // for(auto & sequence : sequences) {
    //     for(auto && window : sequence | ) {
    //         if(window.minimiser_value == max_minimizer) {
    //             kmer = window.window_value ;
    //             break;
    //         }
    //     }
    // }

    return max_occurences;
}


int main(int argc, char** argv)
{
    sharg::parser parser{"ministats", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (sharg::parser_error const &ext) {
        return -1;
    }

    std::vector<seqan3::bitpacked_sequence<seqan3::dna4>> text;
    load_file(args.i, text);
    uint64_t max_minimizer, kmer;
    uint64_t max_occs = count_minimizers(text, args.m, args.k, 1, max_minimizer, kmer);

    std::cout << "minimizer " << decode_kmer(max_minimizer, args.m) << " occurred " << max_occs << " times\n";

    return 0;
}