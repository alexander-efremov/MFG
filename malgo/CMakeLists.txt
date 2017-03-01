cmake_minimum_required(VERSION 2.8)

project(malgo)

add_definitions(-std=c++0x)

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

set(SOURCE_FILES malgo.h)

add_library(malgo STATIC ${SOURCE_FILES} ${HEADER_FILES})
set_target_properties(malgo PROPERTIES LINKER_LANGUAGE C)
