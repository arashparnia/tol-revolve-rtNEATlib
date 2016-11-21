cmake_minimum_required(VERSION 3.6)
project(tol_revolve_rtNEATlib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        main.cpp
        gene.h
        genome.h
        innovation.h
        link.h
        neat.h
        network.h
        nnode.h
        organism.h
        population.h
        species.h
        trait.h
        gene.cpp
        genome.cpp
        innovation.cpp
        link.cpp
        neat.cpp
        network.cpp
        nnode.cpp
        organism.cpp
        population.cpp
        species.cpp
        trait.cpp
        )

add_executable(tol_revolve_rtNEATlib ${SOURCE_FILES})