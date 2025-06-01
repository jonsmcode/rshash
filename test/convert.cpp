#include <iostream>
#include <fstream>
#include <string>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.fasta output.fastq\n";
        return 1;
    }

    std::ifstream fasta(argv[1]);
    std::ofstream fastq(argv[2]);
    if (!fasta.is_open() || !fastq.is_open()) {
        std::cerr << "Error opening files.\n";
        return 1;
    }

    std::string line, header, sequence;
    while (std::getline(fasta, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                // Write previous record
                fastq << '@' << header << '\n'
                      << sequence << '\n'
                      << "+\n"
                      << std::string(sequence.length(), 'I') << '\n';
            }
            header = line.substr(1);
            sequence.clear();
        } else {
            sequence += line;
        }
    }
    // Write last record
    if (!header.empty()) {
        fastq << '@' << header << '\n'
              << sequence << '\n'
              << "+\n"
              << std::string(sequence.length(), 'I') << '\n';
    }

    fasta.close();
    fastq.close();
    return 0;
}