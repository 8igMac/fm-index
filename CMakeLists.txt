cmake_minimum_required(VERSION 3.12)

set(PROJECT_NAME fm-index)
project(${PROJECT_NAME})

# Compiler setup
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 --coverage")

include_directories(${PROJECT_SOURCE_DIR}/include)

# Testing
macro(pkg_add_test TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    target_link_libraries(${TESTNAME} gtest gmock gtest_main)
    add_test(${TESTNAME} ${TESTNAME})
endmacro()
option(BUILD_TESTS "Build tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(${PROJECT_SOURCE_DIR}/submodules/googletest)
    include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    pkg_add_test(saca_k_test unit_test/saca_k_test.cpp)
    pkg_add_test(integration_test unit_test/integration_test.cpp)
endif()

# Regular source file
add_executable(build_sa src/build_sa.cpp)
add_executable(build_index src/build_index.cpp)

# Message
message("Build type: ${CMAKE_BUILD_TYPE}")
message("Build test: ${BUILD_TESTS}")
# message("CMAKE_CXX_FLAGS_DEBUG is ${CMAKE_CXX_FLAGS_DEBUG}")
# message("CMAKE_CXX_FLAGS_RELEASE is ${CMAKE_CXX_FLAGS_RELEASE}")
