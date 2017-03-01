cmake_minimum_required(VERSION 2.8)

project(mfg)

add_definitions(-std=c++0x)

set(HEADER_FILES ${CMAKE_CURRENT_SOURCE_DIR})

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(SOURCE_FILES consts.cpp malgo.h)

add_library(mfg STATIC ${SOURCE_FILES} ${HEADER_FILES})

if(OPENMP_FOUND)
    target_link_libraries(mfg ${OpenMP_CXX_FLAGS})
endif()

