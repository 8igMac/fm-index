cmake_minimum_required(VERSION 3.12)
project(fm-index)
option(BUILD_TESTS "Build tests" ON)

# Dependencies
if(BUILD_TESTS)
    enable_testing()
    include(GoogleTest)
endif()

# Compiler 
set(CMAKE_CXX_STANDARD 14)
include_directories(include)

# Build 
if(BUILD_TESTS)
    add_executable(sais_test unit_test/sais_test.cpp)
    target_link_libraries(sais_test gtest gtest_main)
    add_test(sais_test sais_test)

    # add_executable(bwt_isfm_test unit_test/bwt_isfm_test.cpp)
    # target_link_libraries(bwt_isfm_test gtest gtest_main)
    # add_test(bwt_isfm_test bwt_isfm_test)
endif()

add_executable(build_sa src/build_sa.cpp)
