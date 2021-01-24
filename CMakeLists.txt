cmake_minimum_required(VERSION 3.14)
project(bolotie)

set (CXX11_FOUND "1")
set (CMAKE_PREFIX_PATH "seqan/util/cmake")
set (SEQAN_INCLUDE_PATH "seqan/include")

# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")

find_package (SeqAn REQUIRED)
find_package (OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})

set(CMAKE_CXX_STANDARD 14)

include(FetchContent)

FetchContent_Declare(
        klib
        GIT_REPOSITORY https://github.com/lh3/ksw2.git
)

FetchContent_MakeAvailable(klib)

message(STATUS "klib source dir: ${klib_SOURCE_DIR}")
message(STATUS "klib binary dir: ${klib_BINARY_DIR}")

include_directories(${klib_SOURCE_DIR})

# HTSLIB
include_directories("${CMAKE_SOURCE_DIR}/include/")
link_directories("${CMAKE_SOURCE_DIR}/include/htslib")

include(ExternalProject)
ExternalProject_Add(htslib
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/include/htslib/
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND autoheader COMMAND autoconf COMMAND ./configure --without-libdeflate --disable-libcurl
        BUILD_COMMAND ${MAKE}
        INSTALL_COMMAND ""
        )

add_executable(bolotie
               src/bolotie.cpp
               src/arg_parse.hpp
               src/arg_parse.cpp
               ${klib_SOURCE_DIR}/ksw2.h
               ${klib_SOURCE_DIR}/ksw2_extz2_sse.c
               ${klib_SOURCE_DIR}/ksw2_extz.c
               ${klib_SOURCE_DIR}/ksw2_gg2_sse.c
               ${klib_SOURCE_DIR}/ksw2_gg.c)
add_dependencies(bolotie htslib)
target_link_libraries (bolotie ${SEQAN_LIBRARIES})
target_link_libraries (bolotie pthread -lhts -lz -llzma -lbz2)
