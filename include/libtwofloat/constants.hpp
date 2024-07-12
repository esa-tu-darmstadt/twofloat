#pragma once

/// \file constants.hpp
/// \brief Defines commonly used constants used in all arithmetics.

#include <cmath>
#include <limits>

namespace twofloat {
namespace algorithms {
namespace details {
/// \brief Provides a constexpr pow function, used to calculate some
/// coefficients.
template <typename T>
T constexpr pow(T base, unsigned exponent) {
  return exponent == 0 ? 1 : base * pow(base, exponent - 1);
}
}  // namespace details

/// \brief Defines commonly used constants used in all arithmetics.
template <typename T>
struct constants {
  /// \brief The number of bits in the mantissa of the underlying floating point
  /// type.
  static const constexpr int t = std::numeric_limits<T>::digits;

  /// \brief The maximum value that can be represented by the floating point
  /// type.
  static const constexpr T M = std::numeric_limits<T>::max();

  /// \brief The number of bits that x2 fits into when splitting x into x1+x2 in
  /// algorithms::Split.
  /// \details Calculated as t/2. See Boldo 2006 Algorithm 1 and 2 for details.
  static const constexpr int SplitS = t / 2;

  /// \brief A constant used to split a floating point number into the sum of
  /// two floating point numbers.
  /// \details Calculated as β^(s) + 1. See Boldo 2006 Algorithm 1
  /// and 2 for explanations.
  static const constexpr T SplitC = details::pow(2.0, SplitS) + 1.0;  //- T{};

  /// \brief The threshold used to prevent overflows when splitting a floating
  /// point number into the sum of two numbers.
  /// \details If the absolute value of the number to split is greater than this
  /// threshold, the number is scaled down by SplitScaleFactor before splitting
  /// and the result is scaled back up by SplitScaleFactorInv. Calculated as
  /// M/C. See Boldo 2006 Chapter 3.4 for details.
  static const constexpr T SplitScaleThreshold = M / SplitC;

  /// \brief A factor used to scale down a floating point number before
  /// splitting it to avoid overflows.
  static const constexpr T SplitScaleDownFactor =
      details::pow(0.5, SplitS + 1);  // - T{};

  /// \brief The factor used to scale back up a floating point number after
  /// splitting it.
  static const constexpr T SplitScaleUpFactor = 1 / SplitScaleDownFactor;
};
}  // namespace algorithms
}  // namespace twofloat