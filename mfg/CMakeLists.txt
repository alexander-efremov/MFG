cmake_minimum_required(VERSION 2.8)

project(mfg)

set(CMAKE_CXX_STANDARD 14)

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES consts.cpp solver1.cpp)

add_library(mfg STATIC ${SOURCE_FILES} ${HEADER_FILES})
