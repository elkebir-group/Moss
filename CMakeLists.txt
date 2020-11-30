cmake_minimum_required(VERSION 3.9)

# Project name
project(moss LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

find_package(HTSlib 1.7 REQUIRED)

set(OBJ
    src/moss.cpp
    src/io/bam_io.cpp
    src/io/bam_io.h
    src/io/vcf_io.h
    src/io/vcf_io.cpp
    src/core/types.h
    src/core/types.cpp
    src/core/calling.cpp
    src/core/calling.h
    src/io/loci.h
    src/io/loci.cpp)

add_executable(moss ${OBJ})
target_include_directories(moss PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries(moss ${HTSlib_LIBRARIES})

set(TEST_DIR test)
set(TEST_FILES
    ${TEST_DIR}/types_tests.cpp
    ${TEST_DIR}/catch.hpp
    ${TEST_DIR}/calling_tests.cpp
    ${TEST_DIR}/bam_io_tests.cpp
    src/core/calling.cpp
    src/core/types.cpp
    src/io/vcf_io.cpp)

add_executable(test ${TEST_FILES})
target_include_directories(test PUBLIC ${HTSlib_INCLUDE_DIRS})
target_link_libraries(test ${HTSlib_LIBRARIES})

install(TARGETS moss
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    )
# install(TARGETS test
#     DESTINATION ${CMAKE_INSTALL_PREFIX}
#     )
