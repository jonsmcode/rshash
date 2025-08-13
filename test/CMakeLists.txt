add_executable (datasets datasets.cpp)
target_link_libraries (datasets seqan3::seqan3)

add_executable(verify verify.cpp)
target_link_libraries (verify seqan3::seqan3)

#add_executable(cmp cmp_results.cpp)
#target_link_libraries (cmp seqan3::seqan3)

#add_executable(shrink shrink_queries.cpp)
#target_link_libraries (shrink seqan3::seqan3)

add_executable(minimizer_stats3 minimizer_stats3.cpp)
target_link_libraries (minimizer_stats3 seqan3::seqan3)

add_executable(convert convert.cpp)

add_executable(dna_to_fasta dna_to_fasta.cpp)

add_executable(checkviews check_views.cpp)
target_link_libraries (checkviews seqan3::seqan3)