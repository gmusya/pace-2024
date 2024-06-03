#include "common/assert.h"
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
      best_intersections = (1ull << 48);
      add_to_answer = 0;
      add_min_bound = 0;
      best_matrix = matrix;
    }

    IntersectionMatrix cost_matrix;

    OrderMatrix matrix;
    Graph graph;
    Graph inverse_graph;

    uint64_t current_intersections;
    uint64_t best_intersections;
    uint64_t lower_bound;
    uint64_t add_to_answer;

    uint64_t add_min_bound;

    std::vector<Edge> edge_order;

    OrderMatrix best_matrix;

    UpdatedEdges AddEdge(int from, int to) {
      UpdatedEdges edges;
      edges.emplace_back(from, to);

      for (Vertex first : inverse_graph[from]) {
        for (Vertex last : graph[to]) {
          if (matrix[first][last] == State::Unknown) {
            edges.emplace_back(first, last);
          }
        }
      }
      for (Vertex first : inverse_graph[from]) {
        if (matrix[first][to] == State::Unknown) {
          edges.emplace_back(first, to);
        }
      }
      for (Vertex last : graph[to]) {
        if (matrix[from][last] == State::Unknown) {
          edges.emplace_back(from, last);
        }
      }

      for (const auto& [x, y] : edges) {
        matrix[x][y] = State::Forward;
        matrix[y][x] = State::Backward;
        graph[x].emplace_back(y);
        inverse_graph[y].emplace_back(x);
        current_intersections += cost_matrix[x][y];
        add_min_bound -= std::min(cost_matrix[x][y], cost_matrix[y][x]);
      }

      return edges;
    }

    void UndoEdge(const UpdatedEdges& edges) {
      for (const auto& [x, y] : edges) {
        matrix[x][y] = State::Unknown;
        matrix[y][x] = State::Unknown;
        graph[x].pop_back();
        inverse_graph[y].pop_back();
        current_intersections -= cost_matrix[x][y];
        add_min_bound += std::min(cost_matrix[x][y], cost_matrix[y][x]);
      }
    }
  };

  void Solve(TreeState& tree_state) {
    if (tree_state.current_intersections + tree_state.add_min_bound >=
        tree_state.best_intersections) {
      return;
    }
    ENSURE_OR_THROW(tree_state.current_intersections + tree_state.add_min_bound +
                            tree_state.add_to_answer >=
                    tree_state.lower_bound);
    bool edge_found = false;

    for (const auto& [u, v] : tree_state.edge_order) {
      if (tree_state.matrix[u][v] != State::Unknown) {
        continue;
      }
      edge_found = true;
      std::vector<Edge> order;
      if (tree_state.cost_matrix[u][v] < tree_state.cost_matrix[v][u]) {
        order.emplace_back(u, v);
        order.emplace_back(v, u);
      } else {
        order.emplace_back(v, u);
        order.emplace_back(u, v);
      }
      for (const Edge& edge_to_add : order) {
        if (tree_state.best_intersections + tree_state.add_to_answer == tree_state.lower_bound) {
          return;
        }
        const auto& [from, to] = edge_to_add;
        auto before1 = tree_state.current_intersections;
        auto before2 = tree_state.add_min_bound;

        UpdatedEdges edges = tree_state.AddEdge(from, to);
        Solve(tree_state);
        tree_state.UndoEdge(edges);

        auto after1 = tree_state.current_intersections;
        auto after2 = tree_state.add_min_bound;

        ENSURE_OR_THROW(before1 == after1);
        ENSURE_OR_THROW(before2 == after2);
      }
    }
    if (!edge_found) {
      ENSURE_OR_THROW(tree_state.add_min_bound == 0);

      std::cerr << "tree_state.add_to_answer = " << tree_state.add_to_answer << "\n";
      std::cerr << "best_intersections = "
                << tree_state.current_intersections + tree_state.add_to_answer << "\n";
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
  // ocm::Task task = ocm::Task::FromStream(is);

  std::cerr << "a_size = " << task.a_size << ", b_size = " << task.b_size
            << ", edges_count = " << task.edges.size() << "\n";
  auto graph = ocm::EdgesToGraph(task.edges, task.b_size);
  auto inter_matrix = ocm::BuildIntersectionMatrix(graph);

  uint32_t lower_bound = 0;
  for (ocm::Vertex u = 0; u < graph.size(); ++u) {
    for (ocm::Vertex v = u + 1; v < graph.size(); ++v) {
      lower_bound += std::min(inter_matrix[u][v], inter_matrix[v][u]);
    }
  }
  std::cerr << "lower_bound = " << lower_bound << '\n';

  ocm::tree::TreeState state(inter_matrix);
  state.lower_bound = lower_bound;

  for (ocm::Vertex u = 0; u < state.matrix.size(); ++u) {
    for (ocm::Vertex v = 0; v < state.matrix.size(); ++v) {
      if (u == v) {
        continue;
      }
      if (inter_matrix[u][v] == inter_matrix[v][u]) {
        state.matrix[u][v] = ocm::tree::State::None;
        if (u < v) {
          state.add_to_answer += inter_matrix[u][v];
        }
      }
    }
  }

  for (ocm::Vertex u = 0; u < state.matrix.size(); ++u) {
    for (ocm::Vertex v = 0; v < state.matrix.size(); ++v) {
      if (inter_matrix[u][v] == 0 && inter_matrix[v][u] > 0) {
        if (state.matrix[u][v] == ocm::tree::State::Unknown) {
          state.AddEdge(u, v);
        }
      }
    }
  }

  for (ocm::Vertex u = 0; u < state.matrix.size(); ++u) {
    for (ocm::Vertex v = u + 1; v < state.matrix.size(); ++v) {
      if (state.matrix[u][v] == ocm::tree::State::Unknown) {
        state.edge_order.emplace_back(u, v);
        state.add_min_bound += std::min(state.cost_matrix[u][v], state.cost_matrix[v][u]);
      }
    }
  }

  state.best_matrix = state.matrix;

  std::cerr << "state.add_min_bound = " << state.add_min_bound << std::endl;

  std::sort(state.edge_order.begin(), state.edge_order.end(),
            [&](const auto& lhs, const auto& rhs) {
              return std::max(state.cost_matrix[lhs.first][lhs.second],
                              state.cost_matrix[lhs.second][lhs.first]) -
                             std::min(state.cost_matrix[lhs.first][lhs.second],
                                      state.cost_matrix[lhs.second][lhs.first]) >
                     std::max(state.cost_matrix[rhs.first][rhs.second],
                              state.cost_matrix[rhs.second][rhs.first]) -
                             std::min(state.cost_matrix[rhs.first][rhs.second],
                                      state.cost_matrix[rhs.second][rhs.first]);
            });

  Solve(state);

  std::cerr << "best = " << state.best_intersections + state.add_to_answer << '\n';
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

  auto positions = ocm::SolutionToPositions(result);
  uint64_t inersections_var1 = 0;
  for (ocm::Vertex u = 0; u < state.best_matrix.size(); ++u) {
    for (ocm::Vertex v = 0; v < state.best_matrix.size(); ++v) {
      if (positions[u] < positions[v]) {
        inersections_var1 += state.cost_matrix[u][v];
      }
    }
  }

  std::cerr << "inersections_var1 = " << inersections_var1 << std::endl;
  auto intersections = ocm::CountIntersections(task, positions);
  std::cerr << "intersections = " << intersections << std::endl;
  ENSURE_OR_THROW(inersections_var1 == intersections);
  ENSURE_OR_THROW(intersections == state.best_intersections + state.add_to_answer);
  for (const auto& v : result) {
    std::cout << v + task.a_size + 1 << '\n';
  }
  return 0;
}
