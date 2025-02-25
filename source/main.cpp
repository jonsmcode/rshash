// #include <filesystem>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "dict.h"


struct cmd_arguments
{
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{};
    uint8_t m{};
};

struct my_traits:seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};


void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args)
{
    parser.add_positional_option(args.cmd, "command options: build, query");
    // parser.add_option(args.t, 't', "text", "provide text file", seqan3::input_file_validator{{"fa","fasta"}});
    // parser.add_option(args.q, 'q', "query", "provide query file", seqan3::input_file_validator{{"fa","fasta"}});
    parser.add_option(args.i, 'i', "input", "provide input file");
    parser.add_option(args.q, 'q', "query", "provide query file");
    parser.add_option(args.d, 'd', "dict", "provide dict file");
    parser.add_option(args.k, 'k', "k-mer", "k-mer length");
    parser.add_option(args.m, 'm', "minimiser", "minimiser length");
}


int main(int argc, char** argv)
{
    seqan3::argument_parser parser{"kmerdict", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const &ext) {
        return -1;
    }

    // std::vector<seqan3::dna4> input{"TCATCAGTAGCTACATTACG"_dna4};
    std::vector<seqan3::dna4> input;

    if(args.cmd == "build") {
        if(!parser.is_option_set('i')) {
            std::cout << "provide input file\n";
            return -1;
        }
        if(!parser.is_option_set('d')) {
            std::cout << "provide dict output file\n";
            return -1;
        }
        if(!parser.is_option_set('k') || !parser.is_option_set('m')) {
            std::cout << "specify k and m\n";
            return -1;
        }
        auto reference_stream = seqan3::sequence_file_input<my_traits>{args.i};
        for (auto & record : reference_stream) {
            input = record.sequence();
            break;
        }
        Dictionary dict(args.k, args.m);
        // seqan3::debug_stream << input;
        dict.build(input);
        dict.save(args.d);
    }
    else if(args.cmd == "query") {
        if(!parser.is_option_set('d')) {
            std::cout << "provide dict file\n";
            return -1;
        }
        // if(!parser.is_option_set('q')) {
        //     std::cout << "provide query file\n";
        //     return -1;
        // }
        std::vector<seqan3::dna4> query{"GTAGCTA"_dna4};
        // load query file into memory
        std::vector<uint64_t> positions;

        Dictionary dict;
        dict.load(args.d);
        dict.streaming_query(input, query, positions);

        for (auto pos : positions)
            std::cout << pos << ' ';
        std::cout << '\n';
    }
    else {
        std::cout << "illegal command\n";
        return -1;
    }
 
    return 0;
}

