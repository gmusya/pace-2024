#pragma once

#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace ocm {

using Vertex = uint32_t;
using Edge = std::pair<Vertex, Vertex>;
using Solution = std::vector<Vertex>;
using Position = uint32_t;
using Positions = std::vector<Position>;

struct Task {
  Vertex a_size;
  Vertex b_size;

  // (i, j) where 0 <= i < a_size, 0 <= j < b_size
  std::vector<Edge> edges;

  static Task FromStream(std::istream& is);

  static Task FromFile(const std::string& path);
};

Positions SolutionToPositions(const Solution& solution);

Solution PositionsToSolution(const Positions& positions);

void SaveSolution(const Task& task, const Positions& positions, std::ostream& os);

uint64_t CountIntersections(const Task& task, const Positions& positions);

}// namespace ocm
