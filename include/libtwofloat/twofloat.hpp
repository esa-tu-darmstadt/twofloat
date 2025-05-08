#pragma once

/// \file twofloat.hpp
/// \brief Implements the basic data structure used in all arithmetics.

#include <initializer_list>
#include <libtwofloat/constants.hpp>
#include <limits>
#include <type_traits>
#include <iostream>

/// \brief A namespace for algorithms that operate on pairs of floating point
namespace twofloat {

/// \brief Represents a floating point number as the sum of two
/// floating point numbers.
/// \details The struct is templated to allow it to represent
/// pairs of different floating point types (e.g. float and double).
template <typename T>
struct two {
  // Make sure that T is a floating point type
  static_assert(
      std::is_floating_point<T>::value,
      "twofloat::two<T> can only be instantiated with a floating point "
      "type.");

  using InnerT = T;

  /// \brief The high word of the sum.
  T h;

  /// \brief The low word of the sum.
  T l;

  /// \brief Default constructor
  constexpr two() : h(), l() {}

  /// \brief Constructs an instance from a single floating point number.
  constexpr two(T h) : h(h), l(0) {}

  /// \brief Constructs an instance from two floating point numbers.
  constexpr two(T h, T l) : h(h), l(l) {}

  /// \brief Evaluates the sum to a single floating point number of the
  /// specified type.
  template <typename U = T>
  constexpr U eval() const {
    return static_cast<U>(h) + static_cast<U>(l);
  }

  /// \brief Explicitly convert to a floating point type.
  template <typename DestT>
  explicit constexpr operator DestT() const {
    return eval<DestT>();
  }

  /// \brief Implements the equality comparison operator
  bool constexpr operator==(const two<T> &b) { return h == b.h && l == b.l; }

  template <typename SrcT>
  static constexpr two<T> from(SrcT x) {
    T h = (T)x;
    T l = (T)(x - (SrcT)h);
    return {h, l};
  }

  /// \brief Implements the negation operator
  two<T> constexpr operator-() { return {-h, -l}; }
};

/*********** Equality Comparisons ************/
template <typename T>
inline constexpr bool operator==(const two<T> &a, T b) {
  return (a.h == b && a.l == (T)0.0);
}

template <typename T>
inline constexpr bool operator==(const two<T> &a, const two<T> &b) {
  return (a.h == b.h && a.l == b.l);
}

template <typename T>
inline constexpr bool operator==(T a, const two<T> &b) {
  return (a == b.h && b.l == (T)0.0);
}

/*********** Greater-Than Comparisons ************/
template <typename T>
inline constexpr bool operator>(const two<T> &a, T b) {
  return (a.h > b || (a.h == b && a.l > (T)0.0));
}

template <typename T>
inline constexpr bool operator>(const two<T> &a, const two<T> &b) {
  return (a.h > b.h || (a.h == b.h && a.l > b.l));
}

template <typename T>
inline constexpr bool operator>(T a, const two<T> &b) {
  return (a > b.h || (a == b.h && b.l < (T)0.0));
}

/*********** Less-Than Comparisons ************/
template <typename T>
inline constexpr bool operator<(const two<T> &a, T b) {
  return (a.h < b || (a.h == b && a.l < (T)0.0));
}

template <typename T>
inline constexpr bool operator<(const two<T> &a, const two<T> &b) {
  return (a.h < b.h || (a.h == b.h && a.l < b.l));
}

template <typename T>
inline constexpr bool operator<(T a, const two<T> &b) {
  return (a < b.h || (a == b.h && b.l > (T)0.0));
}

/*********** Greater-Than-Or-Equal-To Comparisons ************/
template <typename T>
inline constexpr bool operator>=(const two<T> &a, T b) {
  return (a.h > b || (a.h == b && a.l >= (T)0.0));
}

template <typename T>
inline constexpr bool operator>=(const two<T> &a, const two<T> &b) {
  return (a.h > b.h || (a.h == b.h && a.l >= b.l));
}

template <typename T>
inline constexpr bool operator>=(T a, const two<T> &b) {
  return (b <= a);
}

/*********** Less-Than-Or-Equal-To Comparisons ************/
template <typename T>
inline constexpr bool operator<=(const two<T> &a, T b) {
  return (a.h < b || (a.h == b && a.l <= (T)0.0));
}

template <typename T>
inline constexpr bool operator<=(const two<T> &a, const two<T> &b) {
  return (a.h < b.h || (a.h == b.h && a.l <= b.l));
}

template <typename T>
inline constexpr bool operator<=(T a, const two<T> &b) {
  return (b >= a);
}

/*********** Not-Equal-To Comparisons ************/
template <typename T>
inline constexpr bool operator!=(const two<T> &a, T b) {
  return (a.h != b || a.l != (T)0.0);
}

template <typename T>
inline constexpr bool operator!=(const two<T> &a, const two<T> &b) {
  return (a.h != b.h || a.l != b.l);
}

template <typename T>
inline constexpr bool operator!=(T a, const two<T> &b) {
  return (a != b.h || b.l != (T)0.0);
}

// Type trait to check if a type is a valid twofloat type
template <typename T>
struct is_twofloat : std::false_type {};

template <typename T>
struct is_twofloat<two<T>> : std::true_type {};

}  // namespace twofloat

namespace std {
/// Absolute value
template <typename T>
inline constexpr twofloat::two<T> abs(twofloat::two<T> a) {
  return (a.h < 0.0) ? -a : a;
}
}  // namespace std