cmake_minimum_required(VERSION 3.10)

project(BurgersEquationSolver VERSION 1.0 DESCRIPTION "Basic 1d Burgers Equation Solver" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable M_PI constant
add_compile_definitions(_USE_MATH_DEFINES)

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/include)

add_compile_options(-Wall -Wextra -Wpedantic)
add_compile_definitions(_USE_MATH_DEFINES)

set(SHARED_SOURCES
    src/burgers.cpp
    src/initial_condition_generator.cpp
    src/schemes/ftcs.cpp
    src/schemes/lax_wendroff.cpp
)

# regular experimental solver
add_executable(burgers
    src/main.cpp
    ${SHARED_SOURCES}
)

# program to create our training data
add_executable(training_data_generator
    src/training_data_generator.cpp
    ${SHARED_SOURCES}
)