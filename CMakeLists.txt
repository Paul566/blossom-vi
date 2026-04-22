cmake_minimum_required(VERSION 3.20)
project(blossom_vi)

include(CheckIPOSupported)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(BLOSSOM_VI_CORE_SOURCES
        MWPMSolver.cpp
        MWPMSolver.h
        DualUpdater.cpp
        DualUpdater.h)

add_executable(blossom_vi main.cpp ${BLOSSOM_VI_CORE_SOURCES})
add_executable(blossom_vi_test
        tests/test_main.cpp
        tests/TesterWeighted.cpp
        tests/TesterWeighted.h
        ${BLOSSOM_VI_CORE_SOURCES})

target_include_directories(blossom_vi_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(blossom_vi PRIVATE $<$<CONFIG:Release>:-march=native>)
    target_compile_options(blossom_vi_test PRIVATE $<$<CONFIG:Release>:-march=native>)
endif()

check_ipo_supported(RESULT ipo_supported)
if(ipo_supported)
    set_property(TARGET blossom_vi PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    set_property(TARGET blossom_vi_test PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
endif()

enable_testing()
add_test(NAME random_correctness COMMAND blossom_vi_test 30 100 20)
