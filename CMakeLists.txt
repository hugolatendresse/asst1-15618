cmake_minimum_required(VERSION 3.10)

project(15618-asmt1)

add_subdirectory(prog1_mandelbrot_threads)
add_subdirectory(prog2_vecintrin)
add_subdirectory(prog3_mandelbrot_ispc)
add_subdirectory(prog4_sqrt)
add_subdirectory(prog5_saxpy)

# Add executable for printpath
# Update the path to printpath.c relative to CMAKE_SOURCE_DIR
message(STATUS "in the root CMAKE_SOURCE_DIR is ${CMAKE_SOURCE_DIR}")
add_executable(printpathroot ${CMAKE_SOURCE_DIR}/printpathroot.c)

## This is a custom target that just builds the printpath executable
#add_custom_target(printpathroot
#        DEPENDS printpathroot
#)
