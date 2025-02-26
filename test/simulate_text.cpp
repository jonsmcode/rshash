#include <filesystem>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>

 using namespace seqan3::literals;


const int k = 20;
const std::vector<seqan3::dna4> repeat{"TCATCAGTAGCTACATTACG"_dna4};

// uint64_t const seed = 0;
uint64_t const seq_length = 100000;
uint64_t const window = 1000;


int main(int argc, char** argv)
{
    std::vector<seqan3::dna4> text;

    int w = window;
    for(int i=0; i < seq_length; i++) {
    	if(w == 0) {
    		// text.insert(i, repeat.begin(), repeat.end());
    		for(int j=0; j < k; j++)
    			text.push_back(repeat[j]);
    		w = window;
    		i+=k;
    	}
    	seqan3::dna4 d;
    	char r = rand() % 4;
    	text.push_back(d.assign_rank(r));
    	w--;
    }

    using types = seqan3::type_list<std::vector<seqan3::dna4>, std::string>;
    using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
    using sequence_record_type = seqan3::sequence_record<types, fields>;

    seqan3::sequence_file_output fout{std::cout, seqan3::format_fasta{}};

    std::string id{"test_id"};
    sequence_record_type record{std::move(text), std::move(id)};

    fout.push_back(record);
    

    return 0;
}

