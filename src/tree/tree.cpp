#include "common/task.h"

#include <algorithm>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>

namespace ocm {

namespace tree {

  enum class State { Unknown, None, Forward, Backward };

  using OrderMatrix = std::vector<std::vector<State>>;
  using UpdatedEdges = std::vector<Edge>;

  struct TreeState {
    TreeState(IntersectionMatrix tmp) : cost_matrix(std::move(tmp)) {
      size_t sz = cost_matrix.size();
      matrix.resize(sz, std::vector<State>(sz));
      graph.resize(sz);
      inverse_graph.resize(sz);
      current_intersections = 0;
      best_intersections = std::numeric_limits<uint64_t>::max();
      best_matrix = matrix;
    }

    IntersectionMatrix cost_matrix;

    OrderMatrix matrix;
    Graph graph;
    Graph inverse_graph;

    uint64_t current_intersections;
    uint64_t best_intersections;
    OrderMatrix best_matrix;
  };

  void Solve(TreeState& tree_state) {
    if (tree_state.current_intersections >= tree_state.best_intersections) {
      return;
    }
    bool edge_found = false;
    for (Vertex u = 0; u < tree_state.matrix.size() && !edge_found; ++u) {
      for (Vertex v = 0; v < tree_state.matrix.size() && !edge_found; ++v) {
        if (tree_state.matrix[u][v] != State::Unknown) {
          continue;
        }
        edge_found = true;
        for (const Edge& edge_to_add : std::vector<Edge>{{u, v}, {v, u}}) {
          const auto& [from, to] = edge_to_add;
          UpdatedEdges edges;
          edges.emplace_back(from, to);
          for (Vertex first : tree_state.inverse_graph[from]) {
            for (Vertex last : tree_state.graph[to]) {
              if (tree_state.matrix[first][last] == State::Unknown) {
                edges.emplace_back(first, last);
              }
            }
          }

          for (const auto& [x, y] : edges) {
            tree_state.matrix[x][y] = State::Forward;
            tree_state.matrix[y][x] = State::Backward;
            tree_state.graph[x].emplace_back(y);
            tree_state.inverse_graph[y].emplace_back(x);
            tree_state.current_intersections += tree_state.cost_matrix[x][y];
          }

          Solve(tree_state);

          for (const auto& [x, y] : edges) {
            tree_state.matrix[x][y] = State::Unknown;
            tree_state.matrix[y][x] = State::Unknown;
            tree_state.graph[x].pop_back();
            tree_state.inverse_graph[y].pop_back();
            tree_state.current_intersections -= tree_state.cost_matrix[x][y];
          }
        }
      }
    }
    if (!edge_found) {
      tree_state.best_matrix = tree_state.matrix;
      tree_state.best_intersections = tree_state.current_intersections;
    }
  }

}// namespace tree

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

  ocm::Task task = ocm::Task::FromStreamCutwidth(is);
  //   ocm::Task task = ocm::Task::FromStream(is);

  std::cerr << "a_size = " << task.a_size << ", b_size = " << task.b_size
            << ", edges_count = " << task.edges.size() << "\n";
  auto graph = ocm::EdgesToGraph(task.edges);
  auto inter_matrix = ocm::BuildIntersectionMatrix(graph);

  uint32_t lower_bound = 0;
  for (ocm::Vertex u = 0; u < graph.size(); ++u) {
    for (ocm::Vertex v = u + 1; v < graph.size(); ++v) {
      lower_bound += std::min(inter_matrix[u][v], inter_matrix[v][u]);
    }
  }
  std::cerr << "lower_bound = " << lower_bound << '\n';

  ocm::tree::TreeState state(inter_matrix);
  Solve(state);

  std::cerr << "best = " << state.best_intersections << '\n';
  std::cerr << "state.matrix.size() = " << state.matrix.size() << '\n';

  std::vector<uint64_t> count(state.matrix.size());
  for (ocm::Vertex u = 0; u < state.best_matrix.size(); ++u) {
    for (ocm::Vertex v = 0; v < state.best_matrix.size(); ++v) {
      if (state.best_matrix[u][v] == ocm::tree::State::Forward) {
        ++count[u];
      }
    }
  }

  std::vector<ocm::Vertex> result(state.matrix.size());
  std::iota(result.begin(), result.end(), 0);

  std::sort(result.begin(), result.end(), [&](const auto& lhs, const auto& rhs) {
    return count[lhs] > count[rhs];
  });

  for (const auto& v : result) {
    std::cout << v + task.a_size + 1 << '\n';
  }
  return 0;
}
