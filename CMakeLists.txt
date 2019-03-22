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
    add_executable(fm_index_test unit_test/fm_index_test.cpp)
    target_link_libraries(fm_index_test gtest gtest_main)
    add_test(fm_index_test fm_index_test)
endif()
