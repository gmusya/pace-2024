#include "common/task.h"
#include "common/assert.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sys/types.h>

namespace ocm {

namespace {
  bool IsPositionsValid(const Task& task, const Positions& positions) {
    if (positions.size() != task.b_size) {
      return false;
    }
    Positions pos = positions;
    std::sort(pos.begin(), pos.end());
    for (Position i = 0; i < pos.size(); ++i) {
      if (pos[i] != i) {
        return false;
      }
    }
    return true;
  }
}// namespace

Task Task::FromStreamCutwidth(std::istream& is) {
  Task result;
  std::string p_token;
  is >> p_token;
  ENSURE_OR_THROW(p_token == "p");
  std::string ocr_token;
  is >> ocr_token;
  ENSURE_OR_THROW(ocr_token == "ocr");
  is >> result.a_size;
  is >> result.b_size;
  uint32_t edges_count;
  is >> edges_count;
  uint32_t cutwidth;
  is >> cutwidth;
  result.cutwidth = cutwidth;
  for (uint32_t i = 0; i < result.a_size + result.b_size; ++i) {
    uint32_t x;
    is >> x;
  }
  result.edges.reserve(edges_count);
  for (uint32_t i = 0; i < edges_count; ++i) {
    Vertex u, v;
    is >> u >> v;
    ENSURE_OR_THROW(1 <= u && u <= result.a_size);
    --u;
    ENSURE_OR_THROW(1 + result.a_size <= v && v <= result.a_size + result.b_size);
    v -= result.a_size + 1;
    result.edges.emplace_back(u, v);
  }
  return result;
}

Task Task::FromStream(std::istream& is) {
  Task result;
  std::string p_token;
  is >> p_token;
  ENSURE_OR_THROW(p_token == "p");
  std::string ocr_token;
  is >> ocr_token;
  ENSURE_OR_THROW(ocr_token == "ocr");
  is >> result.a_size;
  is >> result.b_size;
  uint32_t edges_count;
  is >> edges_count;
  result.edges.reserve(edges_count);
  for (uint32_t i = 0; i < edges_count; ++i) {
    Vertex u, v;
    is >> u >> v;
    ENSURE_OR_THROW(1 <= u && u <= result.a_size);
    --u;
    ENSURE_OR_THROW(1 + result.a_size <= v && v <= result.a_size + result.b_size);
    v -= result.a_size + 1;
    result.edges.emplace_back(u, v);
  }
  return result;
}

Task Task::FromFile(const std::string& path) {
  std::ifstream is(path);
  return FromStream(is);
}

void SaveSolution(const Task& task, const Positions& positions, std::ostream& os) {
  ENSURE_OR_THROW(IsPositionsValid(task, positions));

  Solution solution = PositionsToSolution(positions);

  for (const auto& v : solution) {
    os << v + task.a_size + 1 << '\n';
  }
}

Positions SolutionToPositions(const Solution& solution) {
  Positions result(solution.size());
  for (Position v = 0; v < solution.size(); ++v) {
    result[solution[v]] = v;
  }
  return result;
}

Solution PositionsToSolution(const Positions& positions) {
  Solution result(positions.size());
  for (Vertex v = 0; v < positions.size(); ++v) {
    result[positions[v]] = v;
  }
  return result;
}

#if 0
// O(E^2)
uint64_t CountIntersections(const Task& task, const Solution& solution) {
  ENSURE_OR_THROW(IsPositionsValid(task, solution));
  uint64_t result = 0;
  for (uint32_t ind1 = 0; ind1 < task.edges.size(); ++ind1) {
    for (uint32_t ind2 = ind1 + 1; ind2 < task.edges.size(); ++ind2) {
      auto edge1 = task.edges[ind1];
      auto edge2 = task.edges[ind2];
      if (edge1.first == edge2.first) {
        continue;
      }
      if (edge1.first > edge2.first) {
        std::swap(edge1, edge2);
      }
      // u < x => ((u, v) intersects (x, y) <=> p(v) > p(y))
      result += solution[edge1.second] > solution[edge2.second];
    }
  }
  return result;
}
#endif

namespace {

  class FenwickTree {
public:
    FenwickTree(uint32_t size) : size_(size), values_(size) {
    }

    void Add(uint32_t pos, int32_t val) {
      while (pos < size_) {
        values_[pos] += val;
        pos |= pos + 1;
      }
    }

    int32_t Get(int32_t pos) {
      int32_t result = 0;
      while (pos != (~0)) {
        result += values_[pos];
        pos &= pos + 1;
        --pos;
      }
      return result;
    }

private:
    uint32_t size_;
    std::vector<int32_t> values_;
  };

}// namespace

// O(E log V)
uint64_t CountIntersections(const Task& task, const Positions& positions) {
  ENSURE_OR_THROW(IsPositionsValid(task, positions));

  std::vector<std::vector<Vertex>> neighbors(task.a_size);
  for (const auto& [u, v] : task.edges) {
    neighbors[u].push_back(positions[v]);
  }

  FenwickTree fenwick_tree(task.b_size);
  uint64_t result = 0;
  for (Vertex u = neighbors.size() - 1; u != (~0); --u) {
    for (Vertex v : neighbors[u]) {
      result += fenwick_tree.Get(v - 1);
    }

    for (Vertex v : neighbors[u]) {
      fenwick_tree.Add(v, 1);
    }
  }
  return result;
}

namespace {
  std::pair<uint64_t, uint64_t>
  CountIntersectionsSwappedNaive(const std::vector<ocm::Vertex>& lhs,
                                 const std::vector<ocm::Vertex>& rhs) {
    uint64_t intersections_before_swap = 0;
    uint64_t intersections_after_swap = 0;
    for (const auto& x : lhs) {
      for (const auto& y : rhs) {
        if (x > y) {
          ++intersections_before_swap;
        }
        if (x < y) {
          ++intersections_after_swap;
        }
      }
    }
    return std::make_pair(intersections_before_swap, intersections_after_swap);
  }

  std::pair<uint64_t, uint64_t> CountIntersectionsSwapped(const std::vector<ocm::Vertex>& lhs,
                                                          const std::vector<ocm::Vertex>& rhs) {
    uint64_t intersections_before_swap = 0;
    uint64_t intersections_after_swap = 0;
    uint64_t same_elements = 0;
    uint32_t left_iter = 0;
    uint32_t right_iter = 0;
    while (left_iter + right_iter < lhs.size() + rhs.size()) {
      if (left_iter < lhs.size() &&
          (right_iter == rhs.size() || lhs[left_iter] <= rhs[right_iter])) {
        intersections_before_swap += right_iter;
        if (right_iter < rhs.size() && lhs[left_iter] == rhs[right_iter]) {
          ++same_elements;
        }
        ++left_iter;
      } else {
        ++right_iter;
      }
    }
    intersections_after_swap = lhs.size() * rhs.size() - intersections_before_swap - same_elements;
    return std::make_pair(intersections_before_swap, intersections_after_swap);
  }
}// namespace

Graph EdgesToGraph(const std::vector<Edge>& edges, uint32_t vertex_count) {
  Graph result(vertex_count);
  for (const auto& [u, v] : edges) {
    if (v >= result.size()) {
      result.resize(v + 1);
    }
    result[v].emplace_back(u);
  }
  return result;
}

IntersectionMatrix BuildIntersectionMatrix(const Graph& graph) {
  IntersectionMatrix result(graph.size(), std::vector<uint32_t>(graph.size()));
  for (Vertex u = 0; u < graph.size(); ++u) {
    for (Vertex v = u + 1; v < graph.size(); ++v) {
      auto [lhs, rhs] = CountIntersectionsSwapped(graph[u], graph[v]);
      result[u][v] = lhs;
      result[v][u] = rhs;
    }
  }
  return result;
}

}// namespace ocm
