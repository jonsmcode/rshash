#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>


struct cmd_arguments {
    std::filesystem::path i{};
    std::filesystem::path p{};
    std::filesystem::path q{};
    uint8_t k{};
};

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void initialise_argument_parser(seqan3::argument_parser &parser, cmd_arguments &args) {
    parser.add_option(args.i, 'i', "input", "provide input file");
    parser.add_option(args.q, 'q', "query", "provide query file");
    parser.add_option(args.p, 'p', "positions", "provide position file");
    parser.add_option(args.k, 'k', "k-mer", "k-mer length");
}

int load_text(const std::filesystem::path &filepath, std::vector<seqan3::dna4> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        std::ranges::move(record.sequence(), std::back_inserter(output));
    }
    return 0;
}

int load_queries(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &queries) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        queries.push_back(std::move(record.sequence()));
    }
    return 0;
}

int load_positions(const std::filesystem::path &filepath, std::vector<std::vector<uint64_t>> &positions) {
    std::ifstream file(filepath);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream line_stream(line);
        std::vector<uint64_t> line_positions;
        uint64_t number;
        while (line_stream >> number) {
            line_positions.push_back(number);
        }
        positions.push_back(line_positions);
    }
    file.close();
    return 0;
}

int check_arguments(seqan3::argument_parser &parser, cmd_arguments &args) {
    if(!parser.is_option_set('i'))
        throw seqan3::user_input_error("provide input file.");
    if(!parser.is_option_set('q'))
        throw seqan3::user_input_error("provide query file.");
    if(!parser.is_option_set('p'))
        throw seqan3::user_input_error("provide position file.");
    if(!parser.is_option_set('k'))
        throw seqan3::user_input_error("specify k");

    return 0;
}


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


int main(int argc, char** argv)
{
    seqan3::argument_parser parser{"verify", argc, argv};
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
    load_text(args.i, text);
    std::vector<std::vector<uint64_t>> positions;
    load_positions(args.p, positions);
    std::vector<std::vector<seqan3::dna4>> queries;
    load_queries(args.q, queries);

    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{args.k});

    int q = 0;
    for(auto query : queries) {
        if(positions[q].size() == 0) {
            q++;
            continue;
        }
        int h = 0;
        int p = positions[q][0];
        for (auto && hash : query | kmer_view) {
            if(h == p) {
                seqan3::dna4_vector kmer = kmer_to_string(hash, args.k);
                for(int i=0; i < args.k; i++) {
                    if(text[positions[q][p]+i] != kmer[i]) {
                        std::cout << "not correct\n";
                        return -1;
                    }
                }
                p++;
            }
            h++;
        }
        q++;
    }
    std::cout << "correct\n";
 
    return 0;
}

