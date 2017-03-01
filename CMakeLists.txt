cmake_minimum_required(VERSION 3.6)
project(MFG)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(MFG ${SOURCE_FILES})

include_directories(mfgl)

add_subdirectory(mfgl)
add_subdirectory(mfgt)