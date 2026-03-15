add_executable (rshash rshash.cpp rshash3.cpp rshash2.cpp rshash1.cpp)
target_compile_options(rshash PRIVATE -O3 -march=native -Wno-interference-size)
target_include_directories(rshash PRIVATE ${CMAKE_SOURCE_DIR}/sux)
target_compile_definitions(rshash PRIVATE PARSEARCH)
# target_include_directories(rshash PRIVATE ${CMAKE_SOURCE_DIR}/lib/gtl)
include(FetchContent)
FetchContent_Declare(
    gtl
    GIT_REPOSITORY https://github.com/greg7mdp/gtl.git
    GIT_TAG        v1.2.0 # adjust tag/branch/commit as needed
)
FetchContent_MakeAvailable(gtl)
target_link_libraries(rshash PRIVATE seqan3::seqan3 sharg::sharg gtl)
