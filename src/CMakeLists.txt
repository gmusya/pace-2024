cmake_minimum_required(VERSION 3.22)
project(pace-2024-ocm)

set(BUILD_SHARED_LIBS OFF)
set(CMAKE_EXE_LINKER_FLAGS "-static")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE Release)

add_subdirectory(common)
add_subdirectory(heuristic)

