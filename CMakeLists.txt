cmake_minimum_required(VERSION 3.30)
project(blossom_vi)

include(CheckIPOSupported)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(blossom_vi main.cpp
        TesterWeighted.h
        TesterWeighted.cpp
        MWPMSolver.cpp
        MWPMSolver.h
        DualUpdater.cpp
        DualUpdater.h)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(blossom_vi PRIVATE $<$<CONFIG:Release>:-march=native>)
endif()

check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET blossom_vi PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()
