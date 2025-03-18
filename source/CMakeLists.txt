add_executable (main main.cpp dict.cpp)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (main seqan3::seqan3)

add_executable (locate locate.cpp dict_locate.cpp)
target_include_directories(locate PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (locate seqan3::seqan3)

add_executable (lookup lookup.cpp dict_lookup.cpp)
target_include_directories(lookup PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (lookup seqan3::seqan3)