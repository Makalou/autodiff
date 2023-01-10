cmake_minimum_required(VERSION 3.24)
project(autodiff)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab lab.cpp lab.cpp)

add_subdirectory(3rd_party/googletest)

add_subdirectory(tests)



