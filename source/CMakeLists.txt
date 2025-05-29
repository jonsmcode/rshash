#find_library(PROFILER_LIB profiler REQUIRED)

#add_executable (main main.cpp dict.cpp)
#target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/sux)
##target_link_libraries (main PRIVATE seqan3::seqan3 ${PROFILER_LIB})
#target_link_libraries (main PRIVATE seqan3::seqan3)

add_executable (locate locate.cpp dict_locate.cpp)
target_include_directories(locate PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (locate PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (locate PRIVATE seqan3::seqan3)

add_executable (lookup lookup.cpp dict_lookup.cpp)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (lookup PRIVATE seqan3::seqan3)

add_executable (lookup_unitigs lookup_unitigs.cpp dict_lookup_unitigs.cpp)
target_include_directories(lookup_unitigs PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (lookup_unitigs PRIVATE seqan3::seqan3)

add_executable (comp_lookup comp_lookup.cpp dict_comp.cpp)
target_include_directories(comp_lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (comp_lookup PRIVATE seqan3::seqan3)