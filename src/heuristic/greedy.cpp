#include "common/assert.h"
#include "common/task.h"

#include "heuristic/local.h"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>

namespace hungarian {

// https://cp-algorithms.com/graph/hungarian-algorithm.html#implementation-of-the-hungarian-algorithm

// minimization
std::vector<int> Hungarian(const std::vector<std::vector<int64_t>>& weights) {
  ENSURE_OR_THROW(!weights.empty() && weights.size() <= weights[0].size());

  const int64_t kInfinity = [&]() {
    int64_t result = std::numeric_limits<int64_t>::min();
    for (const auto& row : weights) {
      for (const auto& value : row) {
        result = std::max(result, value);
      }
    }
    return result;
  }();

  int n = static_cast<int>(weights.size()) - 1;
  int m = static_cast<int>(weights[0].size()) - 1;
  std::vector<int> u(n + 1), v(m + 1), p(m + 1), way(m + 1);
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::vector<int> minv(m + 1, kInfinity);
    std::vector<bool> used(m + 1, false);
    do {
      used[j0] = true;
      int i0 = p[j0], delta = kInfinity, j1;
      for (int j = 1; j <= m; ++j)
        if (!used[j]) {
          int cur = weights[i0][j] - u[i0] - v[j];
          if (cur < minv[j])
            minv[j] = cur, way[j] = j0;
          if (minv[j] < delta)
            delta = minv[j], j1 = j;
        }
      for (int j = 0; j <= m; ++j)
        if (used[j])
          u[p[j]] += delta, v[j] -= delta;
        else
          minv[j] -= delta;
      j0 = j1;
    } while (p[j0] != 0);
    do {
      int j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0);
  }

  std::vector<int> ans(n);
  for (int j = 1; j <= m; ++j) {
    ans[p[j] - 1] = j - 1;
  }
  int64_t cost = -v[0];
  std::cerr << "hungarian cost = " << cost << '\n';

  return ans;
}

}// namespace hungarian

namespace ocm {

namespace heuristic {

  constexpr uint64_t kMemoryLimitBytes = 8ull * 1024 * 1024 * 1024;

  namespace greedy {
    using Weight = int64_t;
    using Matrix = std::vector<std::vector<Weight>>;

    std::vector<Weight> BuildWeights(const Task& task,
                                     const std::vector<std::vector<Vertex>>& neighbors, Vertex v) {
      std::vector<Weight> result(task.b_size);
      Weight current_weight = std::accumulate(neighbors[v].begin(), neighbors[v].end(), 0);
      Vertex left_from_me = 0;
      for (Position position = 0; position < task.b_size; ++position) {
        result[position] = current_weight;
        if (left_from_me < neighbors[v].size() && position == neighbors[v][left_from_me]) {
          ++left_from_me;
        }
        current_weight += left_from_me;
        current_weight -= neighbors[v].size() - left_from_me;
      }
      return result;
    }

    Matrix BuildMatrix(const Task& task) {
      ENSURE_OR_THROW(static_cast<uint64_t>(task.b_size) * task.b_size * sizeof(Weight) * 2 <
                      kMemoryLimitBytes);
      Matrix matrix(task.b_size, std::vector<Weight>(task.b_size));

      std::vector<std::vector<Vertex>> neighbors(task.b_size);
      for (const auto& [u, v] : task.edges) {
        neighbors[v].emplace_back(u);
      }

      for (Vertex v = 0; v < task.b_size; ++v) {
        std::sort(neighbors[v].begin(), neighbors[v].end());
        matrix[v] = BuildWeights(task, neighbors, v);
      }

      return matrix;
    }

    Positions Solve(const Task& task) {
      if (task.b_size < 6000) {
        auto matrix = BuildMatrix(task);
        std::vector<std::vector<int64_t>> weights(task.b_size + 1,
                                                  std::vector<int64_t>(task.b_size + 1));
        for (Vertex v = 0; v < task.b_size; ++v) {
          for (Position p = 0; p < task.b_size; ++p) {
            weights[v + 1][p + 1] = matrix[v][p];
          }
        }

        auto solution_hungraian = hungarian::Hungarian(weights);
        Positions positions(solution_hungraian.begin(), solution_hungraian.end());
        std::cerr << "before_optimizations = " << CountIntersections(task, positions) << "\n";
        positions = local::Optimize(task, positions);
        std::cerr << "after_optimizations = " << CountIntersections(task, positions) << "\n";

        return positions;
      }

      std::vector<std::vector<Vertex>> neighbors(task.b_size);
      for (const auto& [u, v] : task.edges) {
        neighbors[v].emplace_back(u);
      }

      std::vector<bool> position_is_taken(task.b_size);
      Positions positions(task.b_size);
      for (Vertex v = 0; v < task.b_size; ++v) {
        std::vector<Vertex> permutation(task.b_size);
        std::iota(permutation.begin(), permutation.end(), 0);
        auto weights = BuildWeights(task, neighbors, v);
        std::sort(permutation.begin(), permutation.end(), [&](Vertex lhs, Vertex rhs) {
          return weights[lhs] < weights[rhs];
        });
        for (Position position : permutation) {
          if (position_is_taken[position]) {
            continue;
          }
          position_is_taken[position] = true;
          positions[v] = position;
          break;
        }
      }
      std::cerr << "before_optimizations = " << CountIntersections(task, positions) << "\n";
      positions = local::Optimize(task, positions);
      std::cerr << "after_optimizations = " << CountIntersections(task, positions) << "\n";
      return positions;
    }


  }// namespace greedy

  namespace trivial {
    Positions Solve(const Task& task) {
      Positions positions(task.b_size);
      std::iota(positions.begin(), positions.end(), 0);
      return positions;
    }
  }// namespace trivial


}// namespace heuristic

}// namespace ocm

int main(int argc, char** argv) {
  std::optional<std::ifstream> input;
  if (argc == 3) {
    input.emplace(argv[1]);
  }
  std::optional<std::ofstream> output;
  if (argc == 3) {
    output.emplace(argv[2]);
  }
  std::istream& is = argc == 3 ? input.value() : std::cin;
  std::ostream& os = argc == 3 ? output.value() : std::cout;

  ocm::Task task = ocm::Task::FromStream(is);
  std::cerr << "a_size = " << task.a_size << ", b_size = " << task.b_size
            << ", edges_count = " << task.edges.size() << "\n";

  uint64_t matrix_size = static_cast<uint64_t>(task.b_size) /* not a typo */ * task.b_size;
  uint64_t bytes = matrix_size * sizeof(ocm::heuristic::greedy::Weight);
  std::cerr << "matrix_size = " << matrix_size << ", matrix_size = " << bytes << "\n";

  ocm::Positions positions;
  positions = ocm::heuristic::greedy::Solve(task);

  std::cerr << "final_score = " << ocm::CountIntersections(task, positions) << "\n";
  ocm::SaveSolution(task, positions, os);
  return 0;
}
