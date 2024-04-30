#include "heuristic/local.h"
#include "common/task.h"
#include <cstdint>
#include <vector>

namespace ocm {

namespace heuristic {

  namespace local {

    namespace {
      std::pair<uint64_t, uint64_t> CountIntersectionsSwapped(const std::vector<Vertex>& lhs,
                                                              const std::vector<Vertex>& rhs) {
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
        intersections_after_swap =
                lhs.size() * rhs.size() - intersections_before_swap - same_elements;
        return std::make_pair(intersections_before_swap, intersections_after_swap);
      }
    }// namespace

    Positions Optimize(const Task& task, const Positions& positions) {
      Solution solution = PositionsToSolution(positions);
      std::vector<std::vector<Vertex>> graph(task.b_size);
      for (const auto& [u, v] : task.edges) {
        graph[v].push_back(u);
      }
      for (Vertex v = 0; v < task.b_size; ++v) {
        std::sort(graph[v].begin(), graph[v].end());
      }
      uint64_t steps = 0;
      while (true) {
        bool optimized = false;
        Position p = 0;
        for (p = 0; p + 1 < positions.size(); ++p) {
          auto pos_lhs = (p + steps) % (positions.size() - 1);
          auto pos_rhs = pos_lhs + 1;
          Vertex v_left = solution[pos_lhs];
          Vertex v_right = solution[pos_rhs];
          auto [before, after] = CountIntersectionsSwapped(graph[v_left], graph[v_right]);
          if (after < before) {
            // auto before_full = CountIntersections(task, SolutionToPositions(solution));
            std::swap(solution[pos_lhs], solution[pos_rhs]);
            optimized = true;
            // auto after_full = CountIntersections(task, SolutionToPositions(solution));
            // ENSURE_OR_THROW(before_full - after_full == before - after);
            break;
          }
        }
        steps += p;
        if (!optimized) {
          break;
        }
      }
      return SolutionToPositions(solution);
    }

  }// namespace local

}// namespace heuristic

}// namespace ocm