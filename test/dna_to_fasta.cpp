// dna_to_fasta.cpp
#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.txt output.fasta\n";
        return 1;
    }

    std::ifstream infile(argv[1]);
    std::ofstream outfile(argv[2]);
    std::string line;
    int seq_num = 1;

    if (!infile) {
        std::cerr << "Error opening input file.\n";
        return 1;
    }
    if (!outfile) {
        std::cerr << "Error opening output file.\n";
        return 1;
    }

    while (std::getline(infile, line)) {
        outfile << ">seq" << seq_num++ << "\n";
        outfile << line << "\n";
    }

    infile.close();
    outfile.close();
    return 0;
}