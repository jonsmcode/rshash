#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "dict.h"


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{};
    uint8_t m{};
};

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args) {
    parser.add_positional_option(args.cmd, "command options: build, query");
    parser.add_option(args.i, 'i', "input", "provide input file");
    parser.add_option(args.q, 'q', "query", "provide query file");
    parser.add_option(args.d, 'd', "dict", "provide dict file");
    parser.add_option(args.k, 'k', "k-mer", "k-mer length");
    parser.add_option(args.m, 'm', "minimiser", "minimiser length");
}

int load_files(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        output.push_back(std::move(record.sequence()));
    }
    return 0;
}

int load_file(const std::filesystem::path &filepath, std::vector<seqan3::dna4> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        std::ranges::move(record.sequence(), std::back_inserter(output));
    }
    return 0;
}

int check_arguments(seqan3::argument_parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw seqan3::user_input_error("provide input file.");
    if(args.cmd == "bq") {
        if(!parser.is_option_set('q'))
            throw seqan3::user_input_error("provide query file.");
        if(!parser.is_option_set('k'))
            throw seqan3::user_input_error("specify k");
    }
    else if(args.cmd == "build") {
        if(!parser.is_option_set('d'))
            throw seqan3::user_input_error("provide dict output file.");
        if(!parser.is_option_set('k'))
            throw seqan3::user_input_error("specify k");
    }
    else if(args.cmd == "query") {
        if(!parser.is_option_set('d'))
            throw seqan3::user_input_error("provide dict file.");
        if(!parser.is_option_set('q'))
            throw seqan3::user_input_error("provide query file.");
    }
    else
        throw seqan3::user_input_error("illegal command");

    return 0;
}


int main(int argc, char** argv)
{
    seqan3::argument_parser parser{"kmerdict", argc, argv};
    cmd_arguments args{};
    initialise_argument_parser(parser, args);
    try {
        parser.parse();
        check_arguments(parser, args);
    }
    catch (seqan3::argument_parser_error const &ext) {
        return -1;
    }

    std::vector<seqan3::dna4> text;
    load_file(args.i, text);

    // if(!parser.is_option_set('m'))
    //     m = ceil(log_4(N)) + 2;
    // else 
    //     m = args.m;

    if(args.cmd == "build") {
        std::cout << "building dict...\n";
        LookupDictionary dict(args.k, args.m);
        dict.build(text);
        std::cout << "done.\n";
        dict.save(args.d);
    }
    else if(args.cmd == "query") {
        LookupDictionary dict;
        dict.load(args.d);

        std::vector<std::vector<seqan3::dna4>> queries;
        load_files(args.q, queries);

        uint64_t n = 0;
        // todo: parallelize queries
        for(auto query : queries) {
            int occurences = dict.streaming_query(text, query);
            // std::cout << occurences << '\n';
            n += occurences;
        }
        std::cout << "total k-mers found: " << n << '\n';
    }
 
    return 0;
}

