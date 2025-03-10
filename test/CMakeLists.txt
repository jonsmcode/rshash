add_executable (bfkmers bfkmers.cpp)
target_link_libraries (bfkmers seqan3::seqan3)

add_executable (datasets datasets.cpp)
target_link_libraries (datasets seqan3::seqan3)