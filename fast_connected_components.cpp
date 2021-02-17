#include <iostream>
#include <vector>
#include <random>
#include <boost/random/taus88.hpp>

template <typename T> std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
  for (auto const& e: v) o << e << " "; return o << "\n";
}

class disjoint_set_t {
  using node_t = struct { size_t p, r, s; };
private:
  mutable std::vector<node_t> v;
public:
  disjoint_set_t() = default;
  disjoint_set_t(disjoint_set_t const&) = default;
  disjoint_set_t(disjoint_set_t&&) = default;
  disjoint_set_t& operator=(disjoint_set_t const&) = default;
  ~disjoint_set_t() = default;
  disjoint_set_t(std::size_t size): v (size) {
    for (std::size_t i = 0; i < size; ++i)
      v[i] = node_t { .p = i, .r = 0, .s = 1 };
  }
  auto find(std::size_t x) const noexcept {
    auto y = x;
    while (v[y].p != y) y = v[y].p;
    while (v[x].p != x) {
      auto& vertex = v[x];
      x = vertex.p;
      vertex.p = y;
    }
    return y;
  }
  auto unite(std::size_t x, std::size_t y) const noexcept {
    auto& rx = v[find(x)];
    auto& ry = v[find(y)];
    if (rx.p != ry.p) {
      if (rx.r < ry.r) { rx.p = ry.p; ry.s += rx.s; }
      else { ry.p = rx.p; rx.s += ry.s; if (rx.r == ry.r) ++rx.r; }
    }
  }
};

class graph_t {
  using edge_t = struct { std::size_t src, dst; };
private:
  std::vector<std::size_t> vertices;
  std::vector<edge_t> edges;
public:
  graph_t() = default;
  graph_t(graph_t const&) = default;
  graph_t(graph_t&&) = default;
  graph_t& operator=(graph_t const&) = default;
  ~graph_t() = default;
  graph_t(std::size_t vs, std::size_t es): vertices(vs), edges(es) {
    for (std::size_t i = 0; i < vs; ++i)
      vertices[i] = i;
    auto complete = std::vector<edge_t> {};
    for (std::size_t i = 0; i < vs; ++i) {
      for (std::size_t j = 0; j < vs; ++j) {
        if (i != j) complete.emplace_back(edge_t {.src = i, .dst = j});
      }
    }
    std::sample(std::begin(complete),
                std::end(complete),
                std::back_inserter(edges),
                es,
                boost::taus88 {std::random_device{}()});
  }

  friend std::ostream& operator<<(std::ostream& o, const graph_t& graph) {
    for (auto const& [src, dst]: graph.edges) o << "[" << src << ", " << dst << "]\n";
    return o;
  }

  auto naive_components() const noexcept {
    auto components = vertices;
    const auto v = std::size(vertices);
    const auto e = std::size(edges);
    for (std::size_t i = 0; i < v; ++i) {
      for (std::size_t j = 0; j < e; ++j) {
        const auto min = std::min(components[edges[j].src], components[edges[j].dst]);
        components[edges[j].src] = min;
        components[edges[j].dst] = min;
      }
    }
    return components;
  }

  auto disjoint_set_component_algorithm() const noexcept {
    auto components = vertices;
    auto disjoint_set = disjoint_set_t(vertices.size());
    const auto v = vertices.size();
    const auto e = edges.size();
    for (size_t i = 0; i < e; ++i) {
      const auto n1 = disjoint_set.find(edges[i].src);
      const auto n2 = disjoint_set.find(edges[i].dst);
      if (n1 != n2)
        disjoint_set.unite(n1, n2);
    }
    for (size_t i = 0; i < v; ++i)
      components[i] = disjoint_set.find(i);
    return components;
  }

};

int main() {
  const auto graph = graph_t(30, 10);
  const auto naive = graph.naive_components();
  const auto fast = graph.disjoint_set_component_algorithm();
  std::cout << naive << std::endl;
  std::cout << fast << std::endl;
  return 0;
}
