#pragma once

#include <cstddef>
#include <functional>
#include <utility>

struct PairHash {
  std::size_t operator()(const std::pair<int, int> &t) const noexcept {
    return std::hash<long long>()(((long long)t.first << 32) ^
                                  (long long)t.second);
  }

  std::size_t operator()(const std::pair<double, double> &t) const noexcept {
    return std::hash<long long>()(((long long)t.first << 32) ^
                                  (long long)t.second);
  }
};
