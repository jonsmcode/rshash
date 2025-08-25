#find_library(PROFILER_LIB profiler REQUIRED)

add_executable (lookup lookup_rsindex.cpp rsindex.cpp rsindex_comp.cpp)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup PRIVATE seqan3::seqan3)

add_executable (lookup2 lookup_rsindex2.cpp rsindex2.cpp rsindex2_comp.cpp)
target_include_directories(lookup2 PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup2 PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup2 PRIVATE seqan3::seqan3 sharg::sharg)

#add_executable (lookup2simd lookup_rsindex2_simd.cpp rsindex2_comp_simd.cpp)
#target_compile_options(lookup3simd PRIVATE -mavx512f -mavx512bw -march=native)
#target_include_directories(lookup3simd PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_include_directories(lookup3simd PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
#target_link_libraries (lookup2simd PRIVATE seqan3::seqan3 sharg::sharg)

add_executable (lookup3 lookup_rsindex3.cpp rsindex3.cpp rsindex3_comp.cpp)
target_include_directories(lookup3 PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_include_directories(lookup3 PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (lookup3 PRIVATE seqan3::seqan3 sharg::sharg)

#add_executable (lookup3simd lookup_rsindex3_simd.cpp rsindex3_comp_simd.cpp)
#target_compile_options(lookup3simd PRIVATE -mavx512f -mavx512bw -march=native)
#target_include_directories(lookup3simd PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_include_directories(lookup3simd PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
#target_link_libraries (lookup3simd PRIVATE seqan3::seqan3 sharg::sharg)

#add_executable (lookup3compsimd lookup_rsindex3_comp_simd.cpp rsindex3_comp3_simd.cpp)
#target_compile_options(lookup3compsimd PRIVATE -mavx512f -mavx512bw -march=native)
#target_include_directories(lookup3compsimd PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_include_directories(lookup3compsimd PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
#target_link_libraries (lookup3compsimd PRIVATE seqan3::seqan3 sharg::sharg)