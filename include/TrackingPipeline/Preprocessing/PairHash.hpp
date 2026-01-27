#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>

struct PairHash {
  std::size_t
  operator()(const std::pair<std::uint16_t, std::uint16_t> &t) const noexcept {
    return std::hash<long long>()(((long long)t.first << 16) ^
                                  (long long)t.second);
  }
};
