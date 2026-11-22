add_executable (hashtable hashtable.cpp)
target_compile_options(hashtable PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(hashtable PRIVATE seqan3::seqan3 sharg::sharg gtl)

# reference for correctness checks: every shape looked up on its own over the query
add_executable (hashtable_ref hashtable_ref.cpp)
target_compile_options(hashtable_ref PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(hashtable_ref PRIVATE seqan3::seqan3 sharg::sharg gtl)

add_executable (datasets test_datasets.cpp)
target_compile_options(datasets PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(datasets PRIVATE seqan3::seqan3 sharg::sharg gtl)

add_executable (ministats mini_stats.cpp)
target_include_directories(ministats PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_compile_definitions(ministats PRIVATE PARSEARCH)
target_compile_options(ministats PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(ministats PRIVATE seqan3::seqan3 sharg::sharg)

add_executable (freqmini freq_mini.cpp)
target_compile_options(freqmini PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(freqmini PRIVATE seqan3::seqan3 sharg::sharg)

add_executable (parse parseinput.cpp)
target_include_directories(parse PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_compile_definitions(parse PRIVATE PARSEARCH)
target_compile_options(parse PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(parse PRIVATE seqan3::seqan3 sharg::sharg gtl)

add_executable (preproc preprocess.cpp)
target_compile_options(preproc PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(preproc PRIVATE seqan3::seqan3 sharg::sharg)
add_executable (samplereads sample_reads.cpp)
target_compile_options(samplereads PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(samplereads PRIVATE seqan3::seqan3 sharg::sharg)

add_executable (eval_locate eval_locate.cpp)
target_compile_options(eval_locate PRIVATE -O3 -march=native)
target_link_libraries(eval_locate PRIVATE sharg::sharg)

# correctness test of the k-mer/shape frequency threshold (rshash build -t) against a reference of all text k-mers
add_executable (test_threshold test_threshold.cpp
    ${CMAKE_SOURCE_DIR}/source/build.cpp ${CMAKE_SOURCE_DIR}/source/lookup.cpp
    ${CMAKE_SOURCE_DIR}/source/locate.cpp ${CMAKE_SOURCE_DIR}/source/io.cpp)
target_include_directories(test_threshold PRIVATE ${CMAKE_SOURCE_DIR}/sux ${CMAKE_SOURCE_DIR}/source)
target_compile_definitions(test_threshold PRIVATE PARSEARCH)
target_compile_options(test_threshold PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(test_threshold PRIVATE seqan3::seqan3 sharg::sharg gtl)

# older generator of artificial datasets (random texts and sampled queries), superseded by datasets (test_datasets.cpp)
add_executable (datasets_random datasets.cpp)
target_compile_options(datasets_random PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(datasets_random PRIVATE seqan3::seqan3)
