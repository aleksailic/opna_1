cmake_minimum_required(VERSION 3.16)
project(opna_1)
set(CMAKE_CXX_STANDARD 17)

include_directories(libs/catch)
include_directories(libs/prettyprint)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

set(FORT_ENABLE_TESTING OFF CACHE INTERNAL "")
add_subdirectory(libs/libfort)

add_executable(opna_1 main.cpp)
target_link_libraries(opna_1 PRIVATE ${Boost_LIBRARIES} fort)

add_executable(opna_1_tests main.cpp)
target_link_libraries(opna_1_tests PRIVATE ${Boost_LIBRARIES} fort)
target_compile_definitions(opna_1_tests PRIVATE TESTING_ENABLED)