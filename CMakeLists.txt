cmake_minimum_required(VERSION 3.30)
project(blossom_vi)

set(CMAKE_CXX_STANDARD 20)

add_executable(blossom_vi main.cpp
        Edge.h
        SolverUnweighted.h
        SolverUnweighted.cpp
        Tester.h
        Tester.cpp
        LabeledDisjointSets.h
        LabeledDisjointSets.cpp
        Node.h
        Node.cpp
        EdgeWeighted.h
        SolverWeighted.h
        SolverWeighted.cpp
        Tree.h
        Tree.cpp)
