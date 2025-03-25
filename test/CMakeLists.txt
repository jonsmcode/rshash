add_executable (datasets datasets.cpp)
target_link_libraries (datasets seqan3::seqan3)

add_executable(verify verify.cpp)
target_link_libraries (verify seqan3::seqan3)