#pragma once

#include <cstddef>
#include <functional>
#include <utility>

namespace detail {

struct PairHash {
  std::size_t operator()(const std::pair<int, int> &t) const noexcept {
    return std::hash<long long>()((static_cast<long long>(t.first) << 32) ^
                                  static_cast<long long>(t.second));
  }

  std::size_t operator()(const std::pair<double, double> &t) const noexcept {
    return std::hash<long long>()((static_cast<long long>(t.first) << 32) ^
                                  static_cast<long long>(t.second));
  }
};

}  // namespace detail
