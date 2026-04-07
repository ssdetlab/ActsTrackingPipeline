#pragma once

#include <vector>

template <typename T>
class ConcatVectorView {
 public:
  ConcatVectorView(const std::vector<T>& a, const std::vector<T>& b)
      : m_a(a), m_b(b) {}

  using value_type = T;

  struct iterator {
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type*;
    using const_pointer = const value_type*;
    using reference = value_type&;
    using const_reference = const value_type&;

    const std::vector<T>* a;
    const std::vector<T>* b;
    difference_type index;

    const_reference operator*() const {
      if (index < a->size()) {
        return (*a)[index];
      } else {
        return (*b)[index - a->size()];
      }
    }

    inline iterator& operator++() noexcept {
      index++;
      return *this;
    }

    inline iterator operator++(int) noexcept {
      iterator tmp = *this;
      index++;
      return tmp;
    }

    inline iterator& operator--() noexcept {
      index--;
      return *this;
    }

    inline iterator operator--(int) noexcept {
      iterator tmp = *this;
      index--;
      return tmp;
    }

    inline iterator& operator+=(difference_type n) noexcept {
      index += n;
      return *this;
    }

    inline iterator& operator-=(difference_type n) noexcept {
      index -= n;
      return *this;
    }

    inline iterator operator+(difference_type n) const noexcept {
      iterator tmp = *this;
      tmp += n;
      return tmp;
    }

    inline iterator operator-(difference_type n) const noexcept {
      iterator tmp = *this;
      tmp -= n;
      return tmp;
    }

    inline difference_type operator-(const iterator& other) const noexcept {
      return index - other.index;
    }

    inline friend iterator operator+(difference_type n,
                                     const iterator& it) noexcept {
      return it + n;
    }

    inline friend iterator operator-(difference_type n,
                                     const iterator& it) noexcept {
      return it + n;
    }

    inline friend bool operator==(const iterator& a, const iterator& b) {
      return (a.index == b.index) && (a.m_typeRange == b.m_typeRange) &&
             (a.m_indexRange == b.m_indexRange);
    };

    inline friend bool operator!=(const iterator& a, const iterator& b) {
      return (a.index != b.index) || (a.m_typeRange != b.m_typeRange) ||
             (a.m_indexRange != b.m_indexRange);
    };
  };

  iterator begin() const { return {&m_a, &m_b, 0}; }

  iterator end() const { return {&m_a, &m_b, m_a.size() + m_b.size()}; }

  iterator::const_reference at(iterator::difference_type i) const {
    return *(begin() + i);
  }

  iterator::difference_type size() const { return m_a.size() + m_b.size(); }

 private:
  const std::vector<T>& m_a;
  const std::vector<T>& m_b;
};
;
