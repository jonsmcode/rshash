#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>

#include "dict4.hpp"


struct cmd_arguments {
    std::string cmd{};
    std::filesystem::path i{};
    std::filesystem::path q{};
    std::filesystem::path d{};
    uint8_t k{};
    uint8_t m{};
};

void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args) {
    parser.add_positional_option(args.cmd, "command options: build, query");
    parser.add_option(args.i, 'i', "input", "provide input file");
    parser.add_option(args.q, 'q', "query", "provide query file");
    parser.add_option(args.d, 'd', "dict", "provide dict file");
    parser.add_option(args.k, 'k', "k-mer", "k-mer length");
    parser.add_option(args.m, 'm', "minimiser", "minimiser length");
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
            throw seqan3::user_input_error("specify m");
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


// void load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
//     auto stream = seqan3::sequence_file_input<my_traits>{filepath};
//     for (auto & record : stream) {
//         output.push_back(std::move(record.sequence()));
//     }
// }

uint64_t load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    uint64_t N = 0;
    for (auto & record : stream) {
        N += record.sequence().size();
        output.push_back(std::move(record.sequence()));
        // if(N >= 1000000000)
        //     break;
        if(N >= 10000000)
            break;
    }
    return N;
}

// todo: ignore queries with 'N's!
// void load_queries(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
//     auto stream = seqan3::sequence_file_input{filepath};
//     for (auto & record : stream) {
//         output.push_back(std::move(record.sequence()));
//     }
// }


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
        UnitigsDictionaryHash2 dict(args.k, args.m);
        dict.build(text);
        std::cout << "done.\n";
        dict.save(args.d);
    }
    else if(args.cmd == "query") {
        std::cout << "loading dict...\n";
        UnitigsDictionaryHash2 dict;
        dict.load(args.d);

        std::cout << "loading queries...\n";
        std::vector<std::vector<seqan3::dna4>> queries;
        load_file(args.q, queries);

        std::cout << "querying...\n";
        uint64_t kmers = 0;
        uint64_t found = 0;
        // todo: parallelize queries
        for(auto query : queries) {
            int occurences = dict.streaming_query(query);
            kmers += query.size()-dict.getk()+1;
            found += occurences;
        }
        std::cout << "==== query report:\n";
        std::cout << "num_kmers = " << kmers << '\n';
        std::cout << "num_positive_kmers = " << found << " (" << (double) found/kmers*100 << "%)\n";
        // int q = 0;
        // for(auto query : queries) {
        //     std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> positions;
        //     dict.streaming_query(query, positions);
        //     for (auto const& [kmer, seq, pos] : positions)
        //         std::cout << '(' << q << ','<< kmer << ',' << seq << ',' << pos << ") ";
        //     q++;
        // }
    }
 
    return 0;
}

