#pragma once

#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>

namespace ocm {

using Vertex = uint32_t;
using Edge = std::pair<Vertex, Vertex>;
using Solution = std::vector<Vertex>;

struct Task {
  Vertex a_size;
  Vertex b_size;

  // (i, j) where 0 <= i < a_size, 0 <= j < b_size
  std::vector<Edge> edges;

  static Task FromStream(std::istream& is);

  static Task FromFile(const std::string& path);
};

uint64_t CountIntersections(const Task& task, const Solution& solution);

}// namespace ocm
