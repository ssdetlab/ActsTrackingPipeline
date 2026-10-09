#pragma once

#include <type_traits>

#include "TrackingPipeline/Utilities/TypeList.hpp"

/// @brief Concept checking whether a type is present in the list
///
/// @tparam T type to check
/// @tparam List list of types to check against
///
/// @return true if T is present in the List, false otherwise
template <typename T, typename List>
concept IsAnyOf = []<typename... Types>(TypeList<Types...>) {
  return (std::is_same_v<Types, T> || ...);
}(List());
