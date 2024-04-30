#include "common/task.h"

#include <fstream>
#include <numeric>
#include <random>

namespace ocm {

namespace heuristic {

  namespace random {
    Solution Solve(const Task& task) {
      Solution solution(task.b_size);
      std::iota(solution.begin(), solution.end(), 0);
      std::mt19937 rng(2101);
      std::shuffle(solution.begin(), solution.end(), rng);
      return solution;
    }
  }// namespace random

}// namespace heuristic

}// namespace ocm

int main(int argc, char** argv) {
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input-file.gr> <output-file.sol>\n";
    return 1;
  }
  std::ifstream is(argv[1]);
  std::ofstream os(argv[2]);
  ocm::Task task = ocm::Task::FromStream(is);
  std::cerr << "a_size = " << task.a_size << ", b_size = " << task.b_size
            << ", edges_count = " << task.edges.size() << "\n";
  std::cerr << "random_score = "
            << ocm::CountIntersections(task, ocm::heuristic::random::Solve(task)) << "\n";
  ocm::SaveSolution(task, ocm::heuristic::random::Solve(task), os);
  return 0;
}