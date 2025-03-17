add_executable (main main.cpp dict.cpp)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_link_libraries (main seqan3::seqan3)
