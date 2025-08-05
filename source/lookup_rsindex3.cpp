#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "rsindex3.hpp"


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{};
    uint8_t m{};
    uint8_t n{};
    uint8_t o{};
};

void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args) {
    parser.add_positional_option(args.cmd, "command options: build, query");
    parser.add_option(args.i, 'i', "input", "provide input file");
    parser.add_option(args.q, 'q', "query", "provide query file");
    parser.add_option(args.d, 'd', "dict", "provide dict file");
    parser.add_option(args.k, 'k', "k-mer", "k-mer length");
    parser.add_option(args.m, 'm', "minimiser1", "minimiser1 length");
    parser.add_option(args.n, 'n', "minimiser2", "minimiser2 length");
    parser.add_option(args.o, 'o', "minimiser3", "minimiser3 length");
}

int check_arguments(seqan3::argument_parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('d'))
            throw seqan3::user_input_error("provide dict file.");
    if(args.cmd == "build") {
        if(!parser.is_option_set('i'))
            throw seqan3::user_input_error("provide input file.");
        if(!parser.is_option_set('k'))
            throw seqan3::user_input_error("specify k");
        if(!parser.is_option_set('m'))
            throw seqan3::user_input_error("specify minimiser1");
        if(!parser.is_option_set('n'))
            throw seqan3::user_input_error("specify minimiser2");
        if(!parser.is_option_set('o'))
            throw seqan3::user_input_error("specify minimiser3");
    }
    else if(args.cmd == "query") {
        if(!parser.is_option_set('q'))
            throw seqan3::user_input_error("provide query file.");
    }
    else
        throw seqan3::user_input_error("illegal command");

    return 0;
}

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    size_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
    }
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

    if(args.cmd == "build") {
        std::cout << "loading text...\n";
        std::vector<std::vector<seqan3::dna4>> text;
        load_file(args.i, text);
        std::cout << "building dict...\n";
        RSIndex dict(args.k, args.m, args.n, args.o);
        dict.build(text);
        std::cout << "done.\n";
        dict.save(args.d);
    }
    else if(args.cmd == "query") {
        std::cout << "loading dict...\n";
        RSIndex dict;
        dict.load(args.d);

        std::cout << "loading queries...\n";
        std::vector<std::vector<seqan3::dna4>> queries;
        load_file(args.q, queries);

        std::cout << "querying...\n";
        uint64_t kmers = 0;
        uint64_t found = 0;
        for(auto query : queries) {
            int occurences = dict.streaming_query(query);
            kmers += query.size()-dict.getk()+1;
            found += occurences;
        }
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
    }
 
    return 0;
}

