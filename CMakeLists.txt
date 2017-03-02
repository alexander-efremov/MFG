cmake_minimum_required(VERSION 3.6)
project(MFG)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp)
add_executable(MFG ${SOURCE_FILES})

include_directories(mfg)
include_directories(malgo)

add_subdirectory(mfg)
add_subdirectory(malgo)
add_subdirectory(tests)

target_link_libraries(MFG mfg)