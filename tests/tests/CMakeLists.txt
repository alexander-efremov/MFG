include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(CMAKE_CXX_STANDARD 14)

set(CPP_SOURCES malgo_tests.cpp mfg_tests.cpp misc_tests.cpp)
add_executable(run_tests ${CPP_SOURCES})

target_link_libraries(run_tests gtest gtest_main mfg)