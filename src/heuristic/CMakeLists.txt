add_executable(heuristic-random main.cpp)
target_link_libraries(heuristic-random common)

add_executable(heuristic-greedy greedy.cpp)
target_link_libraries(heuristic-greedy common)
