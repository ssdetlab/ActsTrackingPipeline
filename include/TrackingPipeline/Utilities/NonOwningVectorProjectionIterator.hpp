#pragma once

#include <cstddef>
#include <iterator>
#include <vector>

template <typename T>
class NonOwningVectorProjectionIterator {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = const value_type*;
  using reference = const value_type&;

  using index_value_type = std::size_t;
  using index_pointer = const index_value_type*;
  using index_reference = const index_value_type&;

  NonOwningVectorProjectionIterator() = delete;

  NonOwningVectorProjectionIterator(const std::vector<T>& typeRange,
                                    const std::vector<std::size_t>& indexRange)
      : m_typeRange(&typeRange),
        m_indexRange(&indexRange),
        m_ptr(indexRange.data()) {};

  reference operator*() const { return m_typeRange->at(*m_ptr); }
  pointer operator->() const { return &m_typeRange->at(*m_ptr); }

  inline NonOwningVectorProjectionIterator& operator++() noexcept {
    m_ptr++;
    return *this;
  }

  inline NonOwningVectorProjectionIterator operator++(int) noexcept {
    NonOwningVectorProjectionIterator tmp = *this;
    m_ptr++;
    return tmp;
  }

  inline NonOwningVectorProjectionIterator& operator--() noexcept {
    m_ptr--;
    return *this;
  }

  inline NonOwningVectorProjectionIterator operator--(int) noexcept {
    NonOwningVectorProjectionIterator tmp = *this;
    m_ptr--;
    return tmp;
  }

  inline NonOwningVectorProjectionIterator& operator+=(
      difference_type n) noexcept {
    m_ptr += n;
    return *this;
  }

  inline NonOwningVectorProjectionIterator& operator-=(
      difference_type n) noexcept {
    m_ptr -= n;
    return *this;
  }

  inline NonOwningVectorProjectionIterator operator+(
      difference_type n) const noexcept {
    NonOwningVectorProjectionIterator tmp = *this;
    tmp += n;
    return tmp;
  }

  inline NonOwningVectorProjectionIterator operator-(
      difference_type n) const noexcept {
    NonOwningVectorProjectionIterator tmp = *this;
    tmp -= n;
    return tmp;
  }

  inline difference_type operator-(
      const NonOwningVectorProjectionIterator& other) const noexcept {
    return m_ptr - other.m_ptr;
  }

  inline friend NonOwningVectorProjectionIterator operator+(
      difference_type n, const NonOwningVectorProjectionIterator& it) noexcept {
    return it + n;
  }

  inline friend NonOwningVectorProjectionIterator operator-(
      difference_type n, const NonOwningVectorProjectionIterator& it) noexcept {
    return it + n;
  }

  inline friend bool operator==(const NonOwningVectorProjectionIterator& a,
                                const NonOwningVectorProjectionIterator& b) {
    return (a.m_ptr == b.m_ptr) && (a.m_typeRange == b.m_typeRange) &&
           (a.m_indexRange == b.m_indexRange);
  };

  inline friend bool operator!=(const NonOwningVectorProjectionIterator& a,
                                const NonOwningVectorProjectionIterator& b) {
    return (a.m_ptr != b.m_ptr) || (a.m_typeRange != b.m_typeRange) ||
           (a.m_indexRange != b.m_indexRange);
  };

 private:
  index_pointer m_ptr;

  const std::vector<T>* m_typeRange;
  const std::vector<std::size_t>* m_indexRange;
};
