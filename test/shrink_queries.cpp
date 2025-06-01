#include <filesystem>
#include <seqan3/io/sequence_file/all.hpp>


void shrink_query(const std::filesystem::path &filepath, const std::string &filename, int n) {
    auto fin = seqan3::sequence_file_input{filepath / filename};
    const std::string fileending = filename.substr(filename.length()-9,filename.length());
    const std::string fileoutname = filename.substr(0,filename.length()-9) + ".K" + std::to_string(n/1000) + fileending;
    auto fout = seqan3::sequence_file_output{filepath / fileoutname};

    int i = 0;
    for(auto & record : fin) {
        if(i == n)
            break;
        fout.push_back(record);
        i++;
    }
}

int main(int argc, char** argv)
{
    if(argc != 4) {
        std::cout << "provide path, filename, n\n";
        return 0;
    }
    const std::filesystem::path path = std::string(argv[1]);
    const std::string filename = std::string(argv[2]);
    int n = atoi(argv[3]);

    shrink_query(path, filename, n);

    return 0;
}

