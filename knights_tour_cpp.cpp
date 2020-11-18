#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <stack>
#include <array>
#include <queue>
#include <functional>
#include <numeric>

constexpr static std::array<std::pair<std::int8_t, std::int8_t>, 8> delta
{{{1, 2}, {1, -2}, {-1, 2}, {-1, -2}, {2, 1}, {2, -1}, {-2, 1}, {-2, -1}}};

auto knights_tour_bfs(
  const std::int32_t knight_x, 
  const std::int32_t knight_y, 
  const std::int32_t target_x,
  const std::int32_t target_y, 
  const std::int32_t max_sz = 8) {

  static const auto isValid = [&](const std::int32_t x, const std::int32_t y) -> bool
  { return (x >= 0 && y >= 0) && (x <= max_sz - 1 && y <= max_sz - 1); };

  if (!isValid(knight_x, knight_y) || !isValid(target_x, target_y))
    return -1;

  std::vector<std::vector<std::int32_t>> grid(max_sz, std::vector<std::int32_t>(max_sz, -1));
  grid[knight_x][knight_y] = 0;

  std::queue<std::int32_t> queue;
  queue.push(knight_x);
  queue.push(knight_y);

  while (!queue.empty()) {
    const auto x = queue.front(); 
    queue.pop();
    const auto y = queue.front();
    queue.pop();
    if ((x == target_x) && (y == target_y)) return grid[x][y];
    for (auto const& d: delta) {
      const auto xn = x + d.first;
      const auto yn = y + d.second;
      if (!isValid(xn, yn)) continue;
      if (grid[xn][yn] == -1) {
        grid[xn][yn] = grid[x][y] + 1;
        queue.push(xn);
        queue.push(yn);
      }
    }
  }
  return -1;
}

auto knights_tour_dfs(
  const std::int32_t knight_x, 
  const std::int32_t knight_y, 
  const std::int32_t target_x, 
  const std::int32_t target_y, 
  const std::int32_t max_sz = 8) {

  std::vector<std::vector<std::int32_t>> grid(max_sz, std::vector<std::int32_t>(max_sz, std::numeric_limits<std::int32_t>::max()));

  static const auto isValid = [&](const std::int32_t x, const std::int32_t y) -> bool
  { return (x >= 0 && y >= 0) && (x <= max_sz - 1 && y <= max_sz - 1); };

  if (!isValid(knight_x, knight_y) || !isValid(target_x, target_y)) 
    return -1;

  std::function<void(const std::int32_t, const std::int32_t, const std::int32_t)> DFS;
  DFS = [&DFS, &max_sz, &grid](const std::int32_t x, const std::int32_t y, const std::int32_t move) {
    if (x < 1 || x > max_sz - 1 || y < 1 || y > max_sz - 1) return;
    if (grid[x][y] <= move)                                 return;
    grid[x][y] = move;
    for (auto const& d: delta) 
      DFS(x - d.first, y - d.second, move + 1);
  };

  DFS(knight_x, knight_y, 0);
  return grid[target_x][target_y];

}

int main() {
  std::cout << "knights_tour_bfs: "               << knights_tour_bfs(5, 7, 1, 1)              << " steps to reach the target.\n";
  std::cout << "knights_tour_dfs: "               << knights_tour_dfs(5, 7, 1, 1)              << " steps to reach the target.\n";
  std::cout << "knights_tour_bfs(10000, 10000): " << knights_tour_bfs(9451, 11, 1, 1, 10000)
                                                                                               << " steps to reach the target.\n";
  std::cout << "knights_tour_dfs(100, 100): "     << knights_tour_dfs(44, 33, 1, 1, 100)       << " steps to reach the target.\n";

  for (;;) {
    std::int32_t kx, ky, tx, ty, sz;
    std::cout << "input: knight_x, knight_y, target_x, target_y, field_size:\n";
    std::cin >> kx >> ky >> tx >> ty >> sz;
    std::cout << "knights_tour_bfs("   << kx  << ", "   << ky << ", "
                                       << tx  << ", "   << ty << ", "
                                       << sz  << "): "  << knights_tour_bfs(kx, ky, tx, ty, sz)
                                       << " steps to reach the target.\n";

    std::cout << "knights_tour_dfs("   << kx  << ", "   << ky << ", "
                                       << tx  << ", "   << ty << ", "
                                       << sz  << "):"   << knights_tour_dfs(kx, ky, tx, ty, sz)
                                       << " steps to reach the target.\n";
  }
  return 0;
}


