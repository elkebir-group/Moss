cmake_minimum_required(VERSION 3.9)

# Project name
project(moss LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/build/cmake/modules/")

find_package(HTSlib 1.4 REQUIRED)

set(OBJ src/main.cpp src/io/bam_io.cpp src/io/bam_io.h src/core/types.h src/core/SnvCaller.cpp src/core/SnvCaller.h src/core/types.cpp)

add_executable(moss ${OBJ})
target_include_directories (moss PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries (moss ${HTSlib_LIBRARIES})