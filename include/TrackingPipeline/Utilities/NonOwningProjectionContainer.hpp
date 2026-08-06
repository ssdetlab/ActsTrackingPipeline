#pragma once

#include <cstddef>
#include <vector>

#include "TrackingPipeline/Utilities/NonOwningProjectionIterator.hpp"

namespace detail {

/// @brief wrapper around the NonOwningProjectionIterator class
/// implementing STL-style access methods
template <typename container_t>
class NonOwningProjectionContainer {
 public:
  using iterator = NonOwningProjectionIterator<container_t>;
  using value_type = iterator::value_type;

  NonOwningProjectionContainer() = delete;

  NonOwningProjectionContainer(const container_t& typeRange,
                               const std::vector<std::size_t>& indexRange)
      : m_typeRange(&typeRange), m_indexRange(&indexRange) {};

  iterator begin() const { return iterator(m_typeRange, m_indexRange); }

  iterator end() const {
    return (iterator(m_typeRange, m_indexRange) + m_indexRange->size());
  }

  iterator::const_reference at(iterator::difference_type i) const {
    return *(begin() + i);
  }

  iterator::difference_type size() const { return m_indexRange->size(); }

 private:
  const container_t* m_typeRange;
  const std::vector<std::size_t>* m_indexRange;
};

}  // namespace detail
