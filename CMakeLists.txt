cmake_minimum_required(VERSION 3.24)
project(autodiff)

set(CMAKE_CXX_STANDARD 20)

add_executable(autodiff main.cpp auto_diff/forward/dual.h auto_diff/auto_diff.h)
