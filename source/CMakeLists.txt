#add_executable (rshash lookup_rsindex3_comp_simd.cpp rsindex3_comp3_simd.cpp rsindex2_comp.cpp rsindex1_comp.cpp)
#target_compile_options(rshash PRIVATE -march=native -Wno-interference-size)
#target_include_directories(rshash PRIVATE ${CMAKE_SOURCE_DIR}/sux)
#target_include_directories(rshash PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
#target_link_libraries (rshash PRIVATE seqan3::seqan3 sharg::sharg)


add_executable (rshashtmp rshash.cpp rshash3.cpp rshash3comp.cpp rshash2.cpp rshash1.cpp)
target_compile_options(rshashtmp PRIVATE -O3 -march=native -Wno-interference-size)
target_include_directories(rshashtmp PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_compile_definitions(rshashtmp PRIVATE PARSEARCH)
target_include_directories(rshashtmp PRIVATE ${CMAKE_SOURCE_DIR}/gtl)
target_link_libraries (rshashtmp PRIVATE seqan3::seqan3 sharg::sharg)
