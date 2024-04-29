add_library(common
    task.cpp
    task.h
)

target_include_directories(common PUBLIC ${CMAKE_SOURCE_DIR})
