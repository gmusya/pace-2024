#include "common/task.h"

#include <numeric>
#include <random>

namespace ocm {

namespace heuristic {

  namespace trivial {
    Solution Solve(const Task& task) {
      Solution solution(task.b_size);
      std::iota(solution.begin(), solution.end(), 0);
      return solution;
    }
  }// namespace trivial

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

int main() {
  ocm::Task task = ocm::Task::FromStream(std::cin);
  std::cerr << "a_size = " << task.a_size << ", b_size = " << task.b_size
            << ", edges_count = " << task.edges.size() << "\n";
  std::cerr << "trivial_score = "
            << ocm::CountIntersections(task, ocm::heuristic::trivial::Solve(task)) << "\n";
  std::cerr << "random_score = "
            << ocm::CountIntersections(task, ocm::heuristic::random::Solve(task)) << "\n";
  ocm::SaveSolution(task, ocm::heuristic::random::Solve(task), std::cout);
  return 0;
}