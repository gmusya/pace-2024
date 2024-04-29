#include "common/task.h"
#include "common/assert.h"

#include <cstdint>
#include <fstream>
#include <sys/types.h>

namespace ocm {

namespace {
  bool IsSolutionValid(const Task& task, const Solution& solution) {
    if (solution.size() != task.b_size) {
      return false;
    }
    Solution sol = solution;
    std::sort(sol.begin(), sol.end());
    for (Vertex i = 0; i < sol.size(); ++i) {
      if (sol[i] != i) {
        return false;
      }
    }
    return true;
  }
}// namespace

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

#if 0
// O(E^2)
uint64_t CountIntersections(const Task& task, const Solution& solution) {
  ENSURE_OR_THROW(IsSolutionValid(task, solution));
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
uint64_t CountIntersections(const Task& task, const Solution& solution) {
  ENSURE_OR_THROW(IsSolutionValid(task, solution));

  std::vector<std::vector<Vertex>> neighbors(task.a_size);
  for (const auto& [u, v] : task.edges) {
    neighbors[u].push_back(solution[v]);
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

}// namespace ocm
