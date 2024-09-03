cmake_minimum_required(VERSION 3.10)

# Project name
project(asmt1-prog2)

# Set the C++ compiler
set(CMAKE_CXX_COMPILER g++)

# Set C++ flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64 -O3 -Wall")

# Include directories
include_directories(../common objs/)

# Directories
set(OBJDIR "objs")
file(MAKE_DIRECTORY ${OBJDIR})

# Source files
message(STATUS "this is the CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
#set(SOURCES
#        main.cpp
#        mandelbrot.cpp
#        ${CMAKE_SOURCE_DIR}/../common/ppm.cpp
#)

## Object files
#set(OBJECTS
#        ${OBJDIR}/main.o
#        ${OBJDIR}/mandelbrot.o
#        ${OBJDIR}/ppm.o
#)
#
## Add executable
#add_executable(${PROJECT_NAME} ${SOURCES})
#
## Link libraries
#target_link_libraries(${PROJECT_NAME} m pthread)
#
## Custom commands for object files
#add_custom_command(
#        OUTPUT ${OBJDIR}/main.o
#        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${CMAKE_SOURCE_DIR}/main.cpp -o ${OBJDIR}/main.o
#        DEPENDS main.cpp ${CMAKE_SOURCE_DIR}/../common/CycleTimer.h
#)
#
#add_custom_command(
#        OUTPUT ${OBJDIR}/mandelbrot.o
#        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${CMAKE_SOURCE_DIR}/mandelbrot.cpp -o ${OBJDIR}/mandelbrot.o
#        DEPENDS mandelbrot.cpp
#)
#
#add_custom_command(
#        OUTPUT ${OBJDIR}/ppm.o
#        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${CMAKE_SOURCE_DIR}/../common/ppm.cpp -o ${OBJDIR}/ppm.o
#        DEPENDS ${CMAKE_SOURCE_DIR}/../common/ppm.cpp
#)

## Clean target
#add_custom_target(clean
#        COMMAND ${CMAKE_COMMAND} -E remove_directory ${OBJDIR}
#        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${PROJECT_NAME}
#        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/*.ppm ${CMAKE_BINARY_DIR}/*~
#)


#
## Add executable for printpath
## Update the path to printpath.c relative to CMAKE_SOURCE_DIR
#add_executable(printpath2 ${CMAKE_SOURCE_DIR}/printpath2.c)
#
## This is a custom target that just builds the printpath executable
#add_custom_target(build_printpath2
#        DEPENDS printpath2
#)