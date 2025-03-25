#find_library(PROFILER_LIB profiler REQUIRED)

add_executable (main main.cpp dict.cpp)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (main PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (main PRIVATE seqan3::seqan3)

add_executable (locate locate.cpp dict_locate.cpp)
target_include_directories(locate PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (locate PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (locate PRIVATE seqan3::seqan3)

add_executable (lookup lookup.cpp dict_lookup.cpp)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (lookup PRIVATE seqan3::seqan3)