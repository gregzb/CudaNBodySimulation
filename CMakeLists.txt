cmake_minimum_required(VERSION 3.17)

project(CudaNBodySimulation VERSION 1.0 LANGUAGES C CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)

set(CMAKE_CXX_FLAGS "-march=native -g -ggdb3")
set(CMAKE_CXX_FLAGS_DEBUG "-march=native -g -ggdb3")

add_compile_options(-lineinfo)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -lineinfo --generate-line-info -lcudadevrt --generate-line-info")

add_subdirectory(src)
add_subdirectory(glad)
add_subdirectory(imgui)