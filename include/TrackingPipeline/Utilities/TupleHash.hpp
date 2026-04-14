#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>

namespace detail {

struct PairHash {
  std::size_t operator()(
      const std::pair<std::uint32_t, std::uint32_t> &t) const noexcept {
    return std::hash<long long>()((static_cast<long long>(t.first) << 32) ^
                                  static_cast<long long>(t.second));
  }

  std::size_t operator()(const std::pair<int, int> &t) const noexcept {
    return std::hash<long long>()((static_cast<long long>(t.first) << 32) ^
                                  static_cast<long long>(t.second));
  }

  std::size_t operator()(const std::pair<double, double> &t) const noexcept {
    return std::hash<long long>()((static_cast<long long>(t.first) << 32) ^
                                  static_cast<long long>(t.second));
  }
};

struct TupleHash {
  std::size_t operator()(
      const std::tuple<std::uint16_t, std::uint16_t, std::uint16_t,
                       std::uint16_t> &t) const noexcept {
    return std::hash<long long>()(
        (static_cast<long long>(std::get<0>(t)) << 48) ^
        (static_cast<long long>(std::get<1>(t)) << 32) ^
        (static_cast<long long>(std::get<2>(t)) << 16) ^
        static_cast<long long>(std::get<3>(t)));
  }
};

}  // namespace detail
