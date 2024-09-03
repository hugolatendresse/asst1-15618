cmake_minimum_required(VERSION 3.10)

# Project name
project(mandelbrot)

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
set(MY_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${MY_SOURCE_DIR}/include)
message(STATUS "prog1 CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "prog1 MY_SOURCE_DIR: ${MY_SOURCE_DIR}")
set(SOURCES
        main.cpp
        mandelbrot.cpp
#        ..common/ppm.cpp
        ${MY_SOURCE_DIR}/../common/ppm.cpp
)

# Object files
set(OBJECTS
        ${OBJDIR}/main.o
        ${OBJDIR}/mandelbrot.o
        ${OBJDIR}/ppm.o
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
target_link_libraries(${PROJECT_NAME} m pthread)

# Custom commands for object files
add_custom_command(
        OUTPUT ${OBJDIR}/main.o
        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${MY_SOURCE_DIR}/main.cpp -o ${OBJDIR}/main.o
        DEPENDS main.cpp ${MY_SOURCE_DIR}/../common/CycleTimer.h
)

add_custom_command(
        OUTPUT ${OBJDIR}/mandelbrot.o
        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${MY_SOURCE_DIR}/mandelbrot.cpp -o ${OBJDIR}/mandelbrot.o
        DEPENDS mandelbrot.cpp
)

add_custom_command(
        OUTPUT ${OBJDIR}/ppm.o
        COMMAND ${CMAKE_CXX_COMPILER} ${CMAKE_CXX_FLAGS} -c ${MY_SOURCE_DIR}/../common/ppm.cpp -o ${OBJDIR}/ppm.o
        DEPENDS ${MY_SOURCE_DIR}/../common/ppm.cpp
)

## Clean target
#add_custom_target(clean
#        COMMAND ${CMAKE_COMMAND} -E remove_directory ${OBJDIR}
#        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/${PROJECT_NAME}
#        COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_BINARY_DIR}/*.ppm ${CMAKE_BINARY_DIR}/*~
#)



# Add executable for printpath
# Update the path to printpath.c relative to MY_SOURCE_DIR
add_executable(printpathprog1 ${MY_SOURCE_DIR}/printpath.c)
