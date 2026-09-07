add_executable (hashtable hashtable.cpp)
target_compile_options(hashtable PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(hashtable PRIVATE seqan3::seqan3 sharg::sharg gtl)

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