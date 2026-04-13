#pragma once

#include <cstddef>
#include <iterator>
#include <vector>

namespace detail {

/// @brief iterator class implementing index projection
/// iteration over a container implementing "at()" method
///
/// @tparam T type of the preimary object vector
///
/// Implements iteration over a container of a primary object T
/// with indices that the iteration is performed over stored in
/// a separate vector.
///
/// Class does not take ownership over the containers with the
/// primary object and the indices, hence lifetime has to be
/// treated carefully
template <typename container_t>
class NonOwningProjectionIterator {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = container_t::value_type;
  using pointer = value_type*;
  using const_pointer = const value_type*;
  using reference = value_type&;
  using const_reference = const value_type&;

  using index_value_type = std::size_t;
  using index_pointer = index_value_type*;
  using const_index_pointer = const index_value_type*;
  using index_reference = index_value_type&;
  using const_index_reference = const index_value_type&;

  NonOwningProjectionIterator() = delete;

  NonOwningProjectionIterator(const container_t& typeRange,
                              const std::vector<std::size_t>& indexRange)
      : m_typeRange(&typeRange),
        m_indexRange(&indexRange),
        m_ptr(indexRange.data()) {};

  NonOwningProjectionIterator(const container_t* typeRange,
                              const std::vector<std::size_t>* indexRange)
      : m_typeRange(typeRange),
        m_indexRange(indexRange),
        m_ptr(indexRange->data()) {};

  const_reference operator*() const { return m_typeRange->at(*m_ptr); }
  const_pointer operator->() const { return &m_typeRange->at(*m_ptr); }

  inline NonOwningProjectionIterator& operator++() noexcept {
    m_ptr++;
    return *this;
  }

  inline NonOwningProjectionIterator operator++(int) noexcept {
    NonOwningProjectionIterator tmp = *this;
    m_ptr++;
    return tmp;
  }

  inline NonOwningProjectionIterator& operator--() noexcept {
    m_ptr--;
    return *this;
  }

  inline NonOwningProjectionIterator operator--(int) noexcept {
    NonOwningProjectionIterator tmp = *this;
    m_ptr--;
    return tmp;
  }

  inline NonOwningProjectionIterator& operator+=(difference_type n) noexcept {
    m_ptr += n;
    return *this;
  }

  inline NonOwningProjectionIterator& operator-=(difference_type n) noexcept {
    m_ptr -= n;
    return *this;
  }

  inline NonOwningProjectionIterator operator+(
      difference_type n) const noexcept {
    NonOwningProjectionIterator tmp = *this;
    tmp += n;
    return tmp;
  }

  inline NonOwningProjectionIterator operator-(
      difference_type n) const noexcept {
    NonOwningProjectionIterator tmp = *this;
    tmp -= n;
    return tmp;
  }

  inline difference_type operator-(
      const NonOwningProjectionIterator& other) const noexcept {
    return m_ptr - other.m_ptr;
  }

  inline friend NonOwningProjectionIterator operator+(
      difference_type n, const NonOwningProjectionIterator& it) noexcept {
    return it + n;
  }

  inline friend NonOwningProjectionIterator operator-(
      difference_type n, const NonOwningProjectionIterator& it) noexcept {
    return it + n;
  }

  inline friend bool operator==(const NonOwningProjectionIterator& a,
                                const NonOwningProjectionIterator& b) {
    return (a.m_ptr == b.m_ptr) && (a.m_typeRange == b.m_typeRange) &&
           (a.m_indexRange == b.m_indexRange);
  };

  inline friend bool operator!=(const NonOwningProjectionIterator& a,
                                const NonOwningProjectionIterator& b) {
    return (a.m_ptr != b.m_ptr) || (a.m_typeRange != b.m_typeRange) ||
           (a.m_indexRange != b.m_indexRange);
  };

 private:
  const_index_pointer m_ptr;

  const container_t* m_typeRange;
  const std::vector<std::size_t>* m_indexRange;
};

}  // namespace detail
