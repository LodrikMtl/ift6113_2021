cmake_minimum_required(VERSION 3.1)
project(hw3)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

set(
        biharmonic_src
        biharmonic.cpp
)

# Add your project files
add_executable(biharmonic ${biharmonic_src})
target_link_libraries(biharmonic igl::core igl::opengl_glfw)

set(
        arap_src
        arap.cpp
)

# Add your project files
add_executable(arap ${arap_src})
target_link_libraries(arap igl::core igl::opengl_glfw)