#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>


const std::filesystem::path textpath = "../../data/salmonella_100_k31_ust.fa.gz";
const std::filesystem::path querypath = "../../data/SRR5833294.10K.fastq.gz";
const int k = 31;
const int m = 15;

uint64_t const seed = 0;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


int load_file(const std::filesystem::path &filepath, std::vector<seqan3::dna4> &output) {
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        output = record.sequence();
        return 0;
    }
    return -1;
}

int load_files(const std::filesystem::path &filepath,
	std::vector<std::vector<seqan3::dna4>> &sequences)
{
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    for (auto & record : stream) {
        sequences.push_back(record.sequence());
    }
    return 0;
}


int main(int argc, char** argv)
{
	std::vector<seqan3::dna4> text;
	std::vector<std::vector<seqan3::dna4>> queries;

	load_file(textpath, text);
	load_files(querypath, queries);

	auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k})
					  | std::views::transform([](uint64_t i) {return i ^ seed;});
	std::vector<uint64_t> positions;

	for(auto query : queries) {
		for (auto && qkmer : query | kmer_view) {
			int i = 0;
			for (auto && tkmer : text | kmer_view) {
				if(qkmer == tkmer)
					positions.push_back(i);
				i++;
			}
		}
	}
	for(auto pos : positions)
        std::cout << pos << ' ';
    std::cout << '\n';
}

