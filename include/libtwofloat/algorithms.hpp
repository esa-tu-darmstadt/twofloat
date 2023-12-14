#pragma once

/// \file algorithms.hpp
/// \brief Implements commonly used algorithms of all arithmetics.

#include <libtwofloat/twofloat.hpp>

namespace twofloat {

/// \brief Provides commonly used algorithms of all algorithms.
namespace algorithms {

/// \brief Fast2Sum algorithm (Dekker 1971).
/// Calculates the sum and resulting rounding error of two floating point
/// numbers using the Fast2Sum algorithm. Requires that |a| >= |b|.
/// Algorithm 1 in https://doi.org/10.1145/3121432
/// \param a The first summand.
/// \param b The second summand.
/// \return The sum of a and b and its error.
template <typename T>
inline two<T> FastTwoSum(T a, T b) {
  two<T> res;
  res.h = a + b;
  T z = res.h - a;
  res.l = b - z;
  return res;
}

/// Calculates the difference and resulting rounding error of two floating point
/// numbers using the same principle as the Fast2Sum algorithm. Requires that
/// |a| >= |b|.
/// \param a The minuend.
/// \param b The subtrahend.
/// \return The difference of a and b and its error.
template <typename T>
inline two<T> FastTwoDiff(T a, T b) {
  two<T> res;
  res.h = a - b;
  T z = a - res.h;
  res.l = z - b;
  return res;
}

/// \brief The TwoSum algorithm (Knuth 1998).
/// Calculates the sum and resulting rounding error of two floating point
/// numbers using the TwoSum algorithm. Algorithm 2 in
/// https://doi.org/10.1145/3121432
/// \param a The first summand.
/// \param b The second summand.
/// \return The sum of a and b and its error .
template <typename T>
inline two<T> TwoSum(T a, T b) {
  two<T> res;
  res.h = a + b;
  T a1 = res.h - b;
  T b1 = res.h - a1;
  T da = a - a1;
  T db = b - b1;
  res.l = da + db;
  return res;
}

/// Calculates the difference and resulting rounding error of two floating point
/// numbers using the same principle as the TwoSum algorithm.
/// \param a The minuend.
/// \param b The subtrahend.
/// \return The difference of a and b and its error .
template <typename T>
inline two<T> TwoDiff(T a, T b) {
  two<T> res;
  res.h = a - b;
  T a1 = res.h - a;
  T b1 = res.h - a1;
  T da = a - b1;
  T db = b + a1;
  res.l = da - db;
  return res;
}

/// \brief Splits a floating point number into the sum of a large and a small
/// floating point number (Veltkamp 1968). Taken from Boldo 2006.
template <typename T>
inline two<T> Split(T x) {
  if (abs(x) > constants<T>::SplitScaleThreshold) {
    // Scale down the number to avoid overflows
    x *= constants<T>::SplitScaleDownFactor;

    T p = x * constants<T>::SplitC;
    T q = x - p;
    T x1 = p + q;
    T x2 = x - x1;

    // Scale the numbers back up
    return {x1 * constants<T>::SplitScaleUpFactor,
            x2 * constants<T>::SplitScaleUpFactor};
  } else {
    // It is safe to split the number without scaling
    T p = x * constants<T>::SplitC;
    T q = x - p;
    T x1 = p + q;
    T x2 = x - x1;
    return {x1, x2};
  }
}

/// \brief Calculates `a*b+c` with infinite precision of the intermediate
/// result.
template <typename T>
inline T fma(T a, T b, T c) {
  if constexpr (std::is_same_v<T, float>)
    return std::fmaf(a, b, c);
  else if constexpr (std::is_same_v<T, double>)
    return std::fma(a, b, c);
  else if constexpr (std::is_same_v<T, long double>)
    return std::fmal(a, b, c);
  else
    static_assert(sizeof(T) == 0, "fma not (yet?) implemented for this type");
}

/// \brief The Fast2Prod algorithm (e.g. Kahan 1996)
/// Calculates the product and resulting rounding error of two floating point
/// numbers using the Fast2Prod algorithm. For performance reasons, this
/// algorithm should only be used if the CPU has support for fused multiply add
/// (fmaf) instructions, which is not the case for most CPUs.
/// Algorithm 3 in https://doi.org/10.1145/3121432
/// \param a The first factor.
/// \param b The second factor.
/// \return The product of a and b and its error.
template <typename T>
inline two<T> Fast2Prod(T a, T b) {
  two<T> res;
  res.h = a * b;
  res.l = fma(a, b, -res.h);
  return res;
}

/// \brief The TwoProd algorithm (Dekker 1971)
/// Calculates the product and resulting rounding error of two floating point
/// numbers using the TwoProd algorithm.
/// \param a The first factor.
/// \param b The second factor.
/// \return The product of a and b and its error.
template <typename T, bool useFMA = false>
inline two<T> TwoProd(T a, T b) {
  if constexpr (useFMA) return Fast2Prod(a, b);

  two<T> res;
  res.h = a * b;

  two<T> a1 = Split(a);
  two<T> b1 = Split(b);

  res.l = ((a1.h * b1.h - res.h) + a1.h * b1.l + a1.l * b1.h) + a1.l * b1.l;
  return res;
}
}  // namespace algorithms
}  // namespace twofloat