cmake_minimum_required(VERSION 3.2)
project(MG)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    addint.c
    interp.c
    main.c
    mg.h
    mglin.c
    nrutil.c
    nrutil.h
    relax.c
    resid.c
    rstrc.c
    slvsml.c)

add_executable(MG ${SOURCE_FILES})