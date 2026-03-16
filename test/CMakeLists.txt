add_executable (hashtable hashtable.cpp)
target_compile_options(hashtable PRIVATE -O3 -march=native -Wno-interference-size)
target_link_libraries(hashtable PRIVATE seqan3::seqan3 sharg::sharg)
