#pragma once

/// \file twofloat.hpp
/// \brief Implements the basic data structure used in all arithmetics.

#include <initializer_list>
#include <libtwofloat/constants.hpp>
#include <limits>
#include <type_traits>

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

  /// \brief The high word of the sum.
  T h;

  /// \brief The low word of the sum.
  T l;

  /// \brief Default constructor
  two() : h(), l() {}

  /// \brief Constructs an instance from a single floating point number.
  two(T h) : h(h), l(0) {}

  /// \brief Constructs an instance from two floating point numbers.
  two(T h, T l) : h(h), l(l) {}

  /// \brief Evaluates the sum to a single floating point number of the
  /// specified type.
  template <typename U = T>
  U eval() const {
    return static_cast<U>(h) + static_cast<U>(l);
  }

  /// \brief Explicitly convert to the underlying floating point type
  explicit operator T() const { return eval(); }

  /// \brief Implements the equality comparison operator
  bool operator==(const two<T> &b) { return h == b.h && l == b.l; }
};

}  // namespace twofloat