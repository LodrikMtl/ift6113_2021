cmake_minimum_required(VERSION 3.1)
project(example)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)

find_package(LIBIGL REQUIRED QUIET)

set(
        main1_src
        main1.cpp
)

# Add your project files
add_executable(main1 ${main1_src})
target_link_libraries(main1 igl::core igl::opengl_glfw)

set(
        main2_src
        main2.cpp
)

# Add your project files
add_executable(main2 ${main2_src})
target_link_libraries(main2 igl::core igl::opengl_glfw)

set(
        main3_src
        main3.cpp
)

# Add your project files
add_executable(main3 ${main3_src})
target_link_libraries(main3 igl::core igl::opengl_glfw)