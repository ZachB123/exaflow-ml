cmake_minimum_required(VERSION 3.14)

project(BurgersEquationSolver VERSION 1.0 DESCRIPTION "Basic 1d Burgers Equation Solver" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable M_PI constant
add_compile_definitions(_USE_MATH_DEFINES)

# nlohmann/json
include(FetchContent)
FetchContent_Declare(json
    URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
    DOWNLOAD_EXTRACT_TIMESTAMP true
)
FetchContent_MakeAvailable(json)

include_directories(${PROJECT_SOURCE_DIR}/include)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

add_compile_definitions(_USE_MATH_DEFINES)

set(SHARED_SOURCES
    src/burgers.cpp
    src/initial_condition_generator.cpp
    src/metadata_writer.cpp
    src/schemes/ftcs.cpp
    src/schemes/lax_wendroff.cpp
    src/schemes/godunov.cpp
)

# regular experimental solver
add_executable(burgers
    src/main.cpp
    ${SHARED_SOURCES}
)
target_link_libraries(burgers PRIVATE nlohmann_json::nlohmann_json)

# program to create our training data
add_executable(training_data_generator
    src/training_data_generator.cpp
    ${SHARED_SOURCES}
)
target_link_libraries(training_data_generator PRIVATE nlohmann_json::nlohmann_json)
