cmake_minimum_required(VERSION 3.22.1)
project(tests)

set(CMAKE_CXX_STANDARD 17)

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

list(APPEND tests
        test_2D_scalar_function.cpp
        test_1D_scalar_function.cpp
        test_1D_scalar_function_optimization.cpp
        test_2D_scalar_function_optimization.cpp
        test_nD_scalar_function.cpp
        test_nD_scalar_function_optimization.cpp)

add_executable(all_tests main.cpp ${tests})

target_link_libraries(all_tests gtest )

include(GoogleTest)
gtest_discover_tests(all_tests)



