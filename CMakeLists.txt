cmake_minimum_required(VERSION 3.30)
project(blossom_vi)

set(CMAKE_CXX_STANDARD 20)

add_executable(blossom_vi main.cpp
        unweighted/Edge.h
        unweighted/SolverUnweighted.h
        unweighted/SolverUnweighted.cpp
        unweighted/Tester.h
        unweighted/Tester.cpp
        unweighted/LabeledDisjointSets.h
        unweighted/LabeledDisjointSets.cpp
        Node.h
        Node.cpp
        EdgeWeighted.h
        EdgeWeighted.cpp
        SolverWeighted.h
        SolverWeighted.cpp
        Tree.h
        Tree.cpp
        TesterWeighted.h
        TesterWeighted.cpp
        Solver.cpp
        Solver.h
        Heap.h
        VzhuhSolver.cpp
        VzhuhSolver.h)
