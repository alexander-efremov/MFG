include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_definitions(-std=c++0x)

set(CPP_SOURCES tests.cpp)
add_executable(run_mfgt ${CPP_SOURCES})

find_package( OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

target_link_libraries(run_mfgt gtest gtest_main)

if(OPENMP_FOUND)
    target_link_libraries(run_mfgt mfg ${OpenMP_CXX_FLAGS})
endif()


