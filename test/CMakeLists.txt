add_executable (bfkmers bfkmers.cpp)
target_link_libraries (bfkmers seqan3::seqan3)

add_executable (simulate_text simulate_text.cpp)
target_link_libraries (simulate_text seqan3::seqan3)

add_executable (simulate_query simulate_query.cpp)
target_link_libraries (simulate_query seqan3::seqan3)