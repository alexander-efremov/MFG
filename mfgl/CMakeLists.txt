cmake_minimum_required(VERSION 2.8)

project(mfg)

add_definitions(-std=c++0x)

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES consts.cpp)

add_library(mfg STATIC ${SOURCE_FILES} ${HEADER_FILES})
