cmake_minimum_required(VERSION 2.8)

#project name
PROJECT(GreenFunction)

message(STATUS "This is BINARY dir " ${CMAKE_BINARY_DIR})
message(STATUS "This is SOURCE dir " ${CMAKE_SOURCE_DIR})

SET(CMAKE_CXX_FLAGS "-std=c++0x")

#head file path
INCLUDE_DIRECTORIES(
)


#library file path
link_directories(
/usr/local/lib/gcc/7
${CMAKE_SOURCE_DIR}/fortran_lib
)


#source directory
add_subdirectory(src)



