#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>


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

int load_file(const std::filesystem::path &filepath, std::vector<std::vector<seqan3::dna4>> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        output.push_back(std::move(record.sequence()));
    }
    return 0;
}

int load_positions(const std::filesystem::path &filepath,
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> &positions)
{
    std::ifstream file(filepath);
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream stream(line);
        std::string tupleStr;

        // Parse each tuple in the line
        while (stream >> tupleStr) {
            // Remove parentheses
            tupleStr.erase(0, 1); // Remove '('
            tupleStr.erase(tupleStr.size() - 1); // Remove ')'

            // Split the tuple into integers
            std::istringstream tupleStream(tupleStr);
            uint64_t a, b, c, d;
            char comma;

            tupleStream >> a >> comma >> b >> comma >> c >> comma >> d;

            // Add the tuple to the positions vector
            positions.emplace_back(a, b, c, d);
        }
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

    std::vector<std::vector<seqan3::dna4>> text;
    load_file(args.i, text);
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t, uint64_t>> positions;
    load_positions(args.p, positions);
    std::vector<std::vector<seqan3::dna4>> queries;
    load_file(args.q, queries);

    
    for(auto const & [q, k, u, i] : positions) {
        for(int j=0; j < args.k; j++) {
            if(queries[q][k + j] != text[u][i+j]) {
                std::cout << "not correct\n";
                return -1;
            }
        }
    }
    
    std::cout << "correct\n";
 
    return 0;
}

