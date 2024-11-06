cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(QrDecompositionSolver VERSION 1.0 LANGUAGES CXX)

add_library(QrLibrary
    lib/argparse.hpp
    lib/csv.hpp
    decomposition.hpp
    matrix.hpp
    reader.hpp
    solver.hpp
)
set_target_properties(QrLibrary PROPERTIES LINKER_LANGUAGE CXX)

add_executable(QrDecompositionSolver main.cpp)
target_link_libraries(QrDecompositionSolver PRIVATE QrLibrary)

if(MSVC)
  target_compile_options(QrDecompositionSolver PRIVATE /W4 /WX)
else()
  target_compile_options(QrDecompositionSolver PRIVATE -Wall -Werror -O3)
endif()
