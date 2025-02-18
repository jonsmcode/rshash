cmake_minimum_required (VERSION 3.20...3.31)
project (seqan3_tutorial CXX)
 
# add seqan3 to search path
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../seqan3/cmake")
 
# require seqan3 with a version between >=3.0.0 and <4.0.0
find_package (seqan3 3.0 REQUIRED)

find_package(Boost 1.40 COMPONENTS program_options REQUIRED)

 
# build app with seqan3
add_executable (dict dict.cpp)
target_link_libraries (dict seqan3::seqan3 Boost::program_options)