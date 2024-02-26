cmake_minimum_required(VERSION 3.14)

project(twofloat VERSION 1.0.0 LANGUAGES CXX)
include(CTest)

add_library(twofloat INTERFACE)
target_include_directories(twofloat INTERFACE include)
target_compile_features(twofloat INTERFACE cxx_std_17)

add_subdirectory(test)