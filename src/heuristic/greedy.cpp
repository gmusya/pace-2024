#include "common/assert.h"
#include "common/task.h"

#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>

namespace ocm {

namespace heuristic {

  constexpr uint64_t kMemoryLimitBytes = 8ull * 1024 * 1024 * 1024;

  namespace greedy {
    using Weight = uint64_t;
    using Matrix = std::vector<std::vector<Weight>>;

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

        Weight current_weight = std::accumulate(neighbors[v].begin(), neighbors[v].end(), 0);
        Vertex left_from_me = 0;
        for (Vertex position = 0; position < task.b_size; ++position) {
          matrix[v][position] = current_weight;
          if (left_from_me < neighbors[v].size() && left_from_me == neighbors[v].size()) {
            ++left_from_me;
          }
          current_weight += left_from_me;
          current_weight -= neighbors[v].size() - left_from_me;
        }
      }

      return matrix;
    }

    Solution Solve(const Task& task) {
      auto matrix = BuildMatrix(task);
      // TODO: use hungarian algorithm

      std::vector<bool> position_is_taken(task.b_size);
      Solution solution(task.b_size);
      for (Vertex v = 0; v < task.b_size; ++v) {
        std::vector<Vertex> permutation(task.b_size);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(), [&](Vertex lhs, Vertex rhs) {
          return matrix[v][lhs] < matrix[v][rhs];
        });
        for (Vertex position : permutation) {
          if (position_is_taken[position]) {
            continue;
          }
          position_is_taken[position] = true;
          solution[v] = position;
          break;
        }
      }
      return solution;
    }


  }// namespace greedy

  namespace trivial {
    Solution Solve(const Task& task) {
      Solution solution(task.b_size);
      std::iota(solution.begin(), solution.end(), 0);
      return solution;
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

  ocm::Solution solution;
  if (bytes * 2 >= ocm::heuristic::kMemoryLimitBytes) {
    std::cerr << "Fallback to trivial solution\n";
    solution = ocm::heuristic::trivial::Solve(task);
  } else {
    solution = ocm::heuristic::greedy::Solve(task);
  }
  std::cerr << "final_score = " << ocm::CountIntersections(task, solution) << "\n";
  ocm::SaveSolution(task, solution, os);
  return 0;
}
