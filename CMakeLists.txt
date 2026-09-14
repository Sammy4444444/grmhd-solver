cmake_minimum_required(VERSION 3.14)
project(grmhd_solver VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization flags for C++ compiler
if(MSVC)
    add_compile_options(/W4 /O2)
else()
    add_compile_options(-Wall -Wextra -O3 -march=native)
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

# Core GRMHD C++ sources
set(SOURCES
    src/main.cpp
)

set(HEADERS
    src/metric.hpp
    src/grid.hpp
    src/hll_solver.hpp
    src/reconstruction.hpp
)

add_executable(grmhd_solver ${SOURCES} ${HEADERS})

message(STATUS "GRMHD Solver configured successfully.")
