cmake_minimum_required (VERSION 3.10)
PROJECT(bvh_distance_queries LANGUAGES CXX CUDA)

include(CheckLanguage)
check_language(CUDA)

FILE(GLOB LIBRARY_FILE "src/*.h" "src/*.hpp" "src/*.cu")

# set(CUDA_ARCHITECTURES 86)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED CONFIG)

include_directories(${pybind11_INCLUDE_DIR})
include_directories(${Python_INCLUDE_DIRS})
include_directories(${Python3_INCLUDE_DIRS})

# if (IS_LINUX)
    # add_library(bvh_distance_queries SHARED ${LIBRARY_FILE})
# else (IS_LINUX)
    add_library(bvh_distance_queries STATIC ${LIBRARY_FILE})
# endif (IS_LINUX)
# set_target_properties(bvh_distance_queries PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

add_executable(bvh_test main.cpp)
target_include_directories(bvh_test PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(bvh_test PUBLIC bvh_distance_queries)
set_target_properties(bvh_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# target_link_libraries(bvh_test PUBLIC CUDA::cudart)
# target_link_libraries(bvh_test PUBLIC CUDA::cufft)

# pybind11_add_module(cuda_distance src/py_binding.cpp)
# target_include_directories(cuda_distance PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
# target_link_libraries(cuda_distance PUBLIC bvh_distance_queries)
# set_target_properties(cuda_distance PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# target_link_libraries(cuda_distance PUBLIC CUDA::cudart)
# target_link_libraries(cuda_distance PUBLIC CUDA::cufft)