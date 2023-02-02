cmake_minimum_required(VERSION 3.22.1)
project(autodiff)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab lab.cpp lab.cpp)

add_subdirectory(3rd_party/googletest)

add_subdirectory(tests)



