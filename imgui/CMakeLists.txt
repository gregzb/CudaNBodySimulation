include_directories(${CudaNBodySimulation_SOURCE_DIR}/include)
file(GLOB imgui_SRC CONFIGURE_DEPENDS "*.h" "*.cpp")
add_library(imgui ${imgui_SRC})
target_link_libraries(imgui PUBLIC glad)