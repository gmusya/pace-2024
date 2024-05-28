#pragma once

#include <cstdint>
#include <iostream>
#include <optional>
#include <utility>
#include <vector>

namespace ocm {

using Vertex = uint32_t;
using Edge = std::pair<Vertex, Vertex>;
using Solution = std::vector<Vertex>;
using Position = uint32_t;
using Positions = std::vector<Position>;
using IntersectionMatrix = std::vector<std::vector<uint32_t>>;
using Graph = std::vector<std::vector<Vertex>>;

struct Task {
  Vertex a_size;
  Vertex b_size;

  // (i, j) where 0 <= i < a_size, 0 <= j < b_size
  std::vector<Edge> edges;

  std::optional<uint32_t> cutwidth;

  static Task FromStream(std::istream& is);

  static Task FromStreamCutwidth(std::istream& is);

  static Task FromFile(const std::string& path);
};

Positions SolutionToPositions(const Solution& solution);

Solution PositionsToSolution(const Positions& positions);

void SaveSolution(const Task& task, const Positions& positions, std::ostream& os);

uint64_t CountIntersections(const Task& task, const Positions& positions);

Graph EdgesToGraph(const std::vector<Edge>& edges);

IntersectionMatrix BuildIntersectionMatrix(const Graph& graph);

}// namespace ocm
