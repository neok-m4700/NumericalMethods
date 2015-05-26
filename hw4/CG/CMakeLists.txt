cmake_minimum_required(VERSION 3.2)
project(CG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    nrutil.h
    nrutil.c
    cg.c)

add_executable(CG ${SOURCE_FILES})