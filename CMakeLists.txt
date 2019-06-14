# cmake_minimum_required(VERSION 3.12)
cmake_minimum_required(VERSION 3.2)
project(fm-index)
option(BUILD_TESTS "Build tests" ON)

# Dependencies
if(BUILD_TESTS)
    enable_testing()
    find_package(GTest REQUIRED)
endif()

# Compiler 
set(CMAKE_CXX_STANDARD 14)
include_directories(include)

# Build 
if(BUILD_TESTS)
    # sais
    add_executable(sais_test unit_test/sais_test.cpp)
    target_link_libraries(sais_test gtest gtest_main)
    add_test(sais_test sais_test)

    # saca-k
    add_executable(saca_k_test unit_test/saca_k_test.cpp)
    target_link_libraries(saca_k_test gtest gtest_main)
    add_test(saca_k_test saca_k_test)

    # # bwt-isfm
    # add_executable(bwt_isfm_test unit_test/bwt_isfm_test.cpp)
    # target_link_libraries(bwt_isfm_test gtest gtest_main)
    # add_test(bwt_isfm_test bwt_isfm_test)
endif()

add_executable(build_sa src/build_sa.cpp)
