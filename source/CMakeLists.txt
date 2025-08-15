#find_library(PROFILER_LIB profiler REQUIRED)


#add_executable (main main.cpp dict.cpp)
#target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/sux)
##target_link_libraries (main PRIVATE seqan3::seqan3 ${PROFILER_LIB})
#target_link_libraries (main PRIVATE seqan3::seqan3)

add_executable (locate locate.cpp dict_locate.cpp)
target_include_directories(locate PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (locate PRIVATE seqan3::seqan3 ${PROFILER_LIB})
target_link_libraries (locate PRIVATE seqan3::seqan3)

#add_executable (lookup lookup.cpp dict_lookup.cpp)
#target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup PRIVATE seqan3::seqan3 ${PROFILER_LIB})
#target_link_libraries (lookup PRIVATE seqan3::seqan3)

add_executable (lookup_unitigs lookup_unitigs.cpp dict_lookup_unitigs.cpp)
target_include_directories(lookup_unitigs PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (lookup_unitigs PRIVATE seqan3::seqan3)

#add_executable (lookup_unitigs lookup_unitigs.cpp dict_unitigs_hash.cpp)
#target_include_directories(lookup_unitigs PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup_unitigs PRIVATE seqan3::seqan3)

#add_executable (lookup_unitigs_simd lookup_unitigs_simd.cpp dict_lookup_unitigs_simd.cpp)
#target_compile_options(lookup_unitigs_simd PRIVATE -mavx512f -mavx512bw -march=native)
#target_include_directories(lookup_unitigs_simd PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup_unitigs_simd PRIVATE seqan3::seqan3)

#add_executable (lookup_unitigsr lookup_unitigsr.cpp dict_unitigs_hash.cpp)
#target_include_directories(lookup_unitigsr PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup_unitigsr PRIVATE seqan3::seqan3)

add_executable (comp_lookup comp_lookup.cpp dict_comp.cpp)
target_include_directories(comp_lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (comp_lookup PRIVATE seqan3::seqan3)

#add_executable (lookup2 lookup_unitigs2.cpp dict_unitigs_hash2.cpp)
#target_include_directories(lookup2 PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_link_libraries (lookup2 PRIVATE seqan3::seqan3)

add_executable (lookup lookup_rsindex.cpp rsindex.cpp rsindex_comp.cpp)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup PRIVATE seqan3::seqan3)

add_executable (lookup2 lookup_rsindex2.cpp rsindex2.cpp rsindex2_comp.cpp)
target_include_directories(lookup2 PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup2 PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup2 PRIVATE seqan3::seqan3)

add_executable (lookup3 lookup_rsindex3.cpp rsindex3.cpp rsindex3_comp.cpp)
target_include_directories(lookup3 PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup3 PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup3 PRIVATE seqan3::seqan3 sharg::sharg)
