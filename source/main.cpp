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
    bool p{false};
};

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
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
    parser.add_flag(args.p, 'p', "positions", "specify if you want to output positions or membership.");
}

int load_file(const std::filesystem::path &filepath, std::vector<seqan3::dna4> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        output = record.sequence();
        return 0;
    }
    return -1;
}

int load_files(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        output.push_back(record.sequence());
    }
    return 0;
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

    std::vector<seqan3::dna4> input;

    if(args.cmd == "bq") {
        if(!parser.is_option_set('i')) {
            std::cout << "provide input file\n";
            return -1;
        }
        if(!parser.is_option_set('q')) {
            std::cout << "provide query file\n";
            return -1;
        }
        if(!parser.is_option_set('k')) {
            std::cout << "specify k\n";
            return -1;
        }
        load_file(args.i, input);
        Dictionary dict(args.k, args.m);
        dict.build(input);

        std::vector<std::vector<seqan3::dna4>> queries;
        load_files(args.q, queries);

        if(args.p) {
            for(auto query : queries) {
                std::vector<uint64_t> positions;
                dict.streaming_query(input, query, positions);
                for (auto pos : positions)
                    std::cout << pos << ' ';
                std::cout << '\n';
            }
        }
        else {
            for(auto query : queries) {
                int occurences = dict.streaming_query(input, query);
                std::cout << occurences << '\n';
            }
        }
    }




    else if(args.cmd == "build") {
        if(!parser.is_option_set('i')) {
            std::cout << "provide input file\n";
            return -1;
        }
        if(!parser.is_option_set('d')) {
            std::cout << "provide dict output file\n";
            return -1;
        }
        if(!parser.is_option_set('k')) {
            std::cout << "specify k\n";
            return -1;
        }
        if(!parser.is_option_set('m')) {
            // m = ceil(log_4(N)) + 2;
        }
        load_file(args.i, input);
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
        if(!parser.is_option_set('q')) {
            std::cout << "provide query file\n";
            return -1;
        }
        Dictionary dict;
        dict.load(args.d);

        std::vector<std::vector<seqan3::dna4>> queries;
        load_files(args.q, queries);
        std::cout << "no queries: " << queries.size() << '\n';

        // todo: parallelize queries
        for(auto query : queries) {
            std::vector<uint64_t> positions;
            dict.streaming_query(input, query, positions);
            for (auto pos : positions)
                std::cout << pos << ' ';
        }
    }
    else {
        std::cout << "illegal command\n";
        return -1;
    }
 
    return 0;
}

