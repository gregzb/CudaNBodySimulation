include_directories(${CudaNBodySimulation_SOURCE_DIR}/include ${CudaNBodySimulation_SOURCE_DIR}/imgui)

find_package(CUDAToolkit)

add_executable(nbodySimulation main.cpp camera.cpp nbody_simulation.cpp nbody_simulation_cuda.cu)
set_target_properties(nbodySimulation PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(nbodySimulation PROPERTIES CUDA_ARCHITECTURES 75)

target_link_libraries(nbodySimulation PRIVATE CUDA::cudart_static
glfw3 GL X11 Xi Xrandr Xxf86vm Xinerama Xcursor rt m pthread dl glad imgui
)

set (source "${CMAKE_SOURCE_DIR}/shaders")
set (destination "${CMAKE_BINARY_DIR}/shaders")
# set (destination "${CMAKE_CURRENT_BINARY_DIR}/shaders")

add_custom_command(
    TARGET nbodySimulation POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${source} ${destination} 
    DEPENDS ${destination}
    COMMENT "symbolic link resources folder from ${source} => ${destination}"
)