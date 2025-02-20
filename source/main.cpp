#include <iostream>
#include "dict.h"


int help(char* arg0) {
    std::cout << "k-mer dictionary"
              << std::endl
              << std::endl;
    std::cout << "Usage: " << arg0 << " <tool> ...\n\n"
              << "Available tools:\n"
              << "  build              \t build a dictionary \n"
              << "  query              \t query a dictionary \n"
              << std::endl;
    return 1;
}


int main(int argc, char** argv)
{
    // if (argc < 2) return help(argv[0]);
    // auto tool = std::string(argv[1]);
    // if (tool == "build") {
    //     // return build(argc-1, argv+1);
    // } else if (tool == "query") {
    //     // return query(argc-1, argv+1);
    // }
    // std::cout << "Unsupported tool '" << tool << "'.\n" << std::endl;
    // return help(argv[0]);

    // auto reference_stream = seqan3::sequence_file_input{reference_file};

    // // read into memory
    // std::vector<std::vector<seqan3::dna4>> text;
    // for (auto& record : reference_stream) {
    //     text.push_back(record.sequence());
    // }
    std::vector<seqan3::dna4> text{"TCATCAGTAGCTACATTACG"_dna4};

    Dictionary dict(5, 2);
    dict.build(text);

    // std::vector<seqan3::dna4> query{"GTAGCTAGCTACA"_dna4};
    std::vector<seqan3::dna4> query{"GTAGCTA"_dna4};
    std::vector<uint64_t> positions;
    dict.streaming_query(text, query, positions);

    for (auto pos : positions)
        std::cout << pos << ' ';
    std::cout << '\n';

    
    return 0;
}