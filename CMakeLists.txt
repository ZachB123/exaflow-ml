cmake_minimum_required(VERSION 3.10)

project(BurgersEquationSolver VERSION 1.0 DESCRIPTION "Basic 1d Burgers Equation Solver" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

add_compile_options(-Wall -Wextra -Wpedantic)

# Source files
set(SOURCES
    src/main.cpp
    src/burgers.cpp
    src/initial_condition_generator.cpp
    src/stencils/ftcs.cpp
    src/stencils/lax_wendroff.cpp
)

# Executable
add_executable(burgers ${SOURCES})