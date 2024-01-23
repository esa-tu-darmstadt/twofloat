#pragma once

/// \file double-word-arithmetic.hpp
/// \brief Implements the double-word arithmetic proposed by Joldeş et al.
/// (2017) for sum and multiplications and Lefèvre et al. (2022) for the square
/// root.

#include <libtwofloat/algorithms.hpp>
#include <libtwofloat/twofloat.hpp>

namespace twofloat {
namespace doubleword {

enum class Mode { Fast, Accurate, Sloppy };

/// \brief Adds a double-word floating point number with a normal floating point
/// number.
/// \details This is algorithm `DWPlusFP` in Joldeş et al. (2017).
template <typename T>
inline two<T> add(const two<T> &x, T y) {
  two<T> s = algorithms::TwoSum(x.h, y);
  T v = x.l + s.l;
  return algorithms::FastTwoSum(s.h, v);
}
template <typename T>
inline two<T> add(T x, const two<T> &y) {
  return add(y, x);
}

/// \brief Subtracts a double-word floating point number with a normal floating
/// point number.
/// \details Derived algorithm `DWPlusFP` in Joldeş et al. (2017).
template <typename T>
inline two<T> sub(const two<T> &x, T y) {
  two<T> s = algorithms::TwoDiff(x.h, y);
  T v = x.l + s.l;
  return algorithms::FastTwoSum(s.h, v);
}
/// \brief Subtracts a double-word floating point number with a normal floating
/// point number.
/// \details Derived algorithm `DWPlusFP` in Joldeş et al. (2017).
template <typename T>
inline two<T> sub(T x, const two<T> &y) {
  two<T> s = algorithms::TwoDiff(x, y.h);
  T v = s.l - y.l;
  return algorithms::FastTwoSum(s.h, v);
}

/// \brief Adds two double-word floating point numbers.
/// \details When using the sloppy mode, the relative error is only bounded if
/// x.h and y.h have the same sign. When using the accurate mode, the relative
/// error is always bounded.
/// \param x The first double-word floating point number.
/// \param y The second double-word floating point number.
/// \tparam mode The mode (sloppy or accurate).
template <Mode mode, typename T>
inline two<T> add(const two<T> &x, const two<T> &y) {
  if constexpr (mode == Mode::Sloppy) {
    // SloppyDWPlusDW in Joldes et al. (2017)
    two<T> s = algorithms::TwoSum(x.h, y.h);
    T v = x.l + y.l;
    T w = s.l + v;
    return algorithms::FastTwoSum(s.h, w);
  } else if constexpr (mode == Mode::Accurate) {
    // AccurateDWPlusDW in Joldes et al. (2017)
    two<T> s = algorithms::TwoSum(x.h, y.h);
    two<T> t = algorithms::TwoSum(x.l, y.l);
    T c = s.l + t.h;
    two<T> v = algorithms::FastTwoSum(s.h, c);
    T w = t.l + v.l;
    return algorithms::FastTwoSum(v.h, w);
  } else
    static_assert(sizeof(T) == 0, "Sloppy and accurate modes are supported");
}

/// \brief Subtracts two double-word floating point numbers.
/// \details When using the sloppy mode, the relative error is only bounded if
/// x.h and y.h have the same sign. When using the accurate mode, the relative
/// error is always bounded.
/// \param x The first double-word floating point number.
/// \param y The second double-word floating point number.
/// \tparam mode The mode (sloppy or accurate).
template <Mode mode, typename T>
inline two<T> sub(const two<T> &x, const two<T> &y) {
  if constexpr (mode == Mode::Sloppy) {
    // Based on SloppyDWPlusDW in Joldes et al. (2017)
    two<T> s = algorithms::TwoDiff(x.h, y.h);
    T v = x.l - y.l;
    T w = s.l + v;
    return algorithms::FastTwoSum(s.h, w);
  } else if constexpr (mode == Mode::Accurate) {
    // AccurateDWPlusDW in Joldes et al. (2017)
    two<T> s = algorithms::TwoDiff(x.h, y.h);
    two<T> t = algorithms::TwoDiff(x.l, y.l);
    T c = s.l + t.h;
    two<T> v = algorithms::FastTwoSum(s.h, c);
    T w = t.l + v.l;
    return algorithms::FastTwoSum(v.h, w);
  } else
    static_assert(sizeof(T) == 0, "Sloppy and accurate modes are supported");
}

/// \brief Multiplies a double-word floating point number with a floating point.
/// \details The accurate algorithm was proposed by Li et al. (2000). The sloppy
/// algorithm was proposed by Higgs (1988).
/// \param x The double-word floating point number.
/// \param y The floating point number.
/// \tparam p The mode (sloppy or accurate). Ignored when using FMA.
/// \tparam useFMA Whether to use FMA instructions.
/// \return The product of x and y.
template <Mode p, bool useFMA, typename T>
inline two<T> mul(const two<T> &x, T y) {
  if constexpr (useFMA) {
    // DWTimesFP3 in Joldes et al. (2017)
    two<T> c = algorithms::Fast2Prod(x.h, y);
    T cl3 = algorithms::fma(x.l, y, c.l);
    return algorithms::FastTwoSum(c.h, cl3);
  }

  if constexpr (p == Mode::Fast) {
    // DWTimesFP2 in Joldes et al. (2017)
    two<T> c = algorithms::TwoProd(x.h, y);
    T cl2 = x.l * y;
    T cl3 = c.l + cl2;
    return algorithms::FastTwoSum(c.h, cl3);
  } else if constexpr (p == Mode::Accurate) {
    // DWTimesFP1 in Joldes et al. (2017)
    two<T> c = algorithms::TwoProd(x.h, y);
    T cl2 = x.l * y;
    two<T> t = algorithms::FastTwoSum(c.h, cl2);
    T tl2 = t.l + c.l;
    return algorithms::FastTwoSum(t.h, tl2);
  } else
    static_assert(sizeof(T) == 0, "Fast, accurate and FMA modes are supported");
}
template <Mode p, bool useFMA, typename T>
inline two<T> mul(T x, const two<T> &y) {
  return mul<p, useFMA>(y, x);
}

/// \brief Multiplies two double-word floating point numbers.
/// \details The non-FMA fast algorithm was proposed by Dekker (1971). The
/// FMA algorithms were proposed by Joldeş et al. (2017).
/// \param x The first double-word floating point number.
/// \param y The second double-word floating point number.
/// \tparam p The mode (fast or accurate). When not using FMA, only the fast
/// mode is supported.
/// \tparam useFMA Whether to use FMA instructions.
template <Mode p, bool useFMA, typename T>
inline two<T> mul(const two<T> &x, const two<T> &y) {
  if constexpr (useFMA) {
    if constexpr (p == Mode::Fast) {
      // DWTimesDW2 in Joldes et al. (2017)
      two<T> c = algorithms::Fast2Prod(x.h, y.h);
      T tl = x.h * y.l;
      T cl2 = algorithms::fma(x.l, y.h, tl);
      T cl3 = c.l + cl2;
      return algorithms::FastTwoSum(c.h, cl3);
    } else if constexpr (p == Mode::Accurate) {
      // DWTimesDW3 in Joldes et al. (2017)
      two<T> c = algorithms::Fast2Prod(x.h, y.h);
      T tl0 = x.l * y.l;
      T tl1 = algorithms::fma(x.h, y.l, tl0);
      T cl2 = algorithms::fma(x.l, y.h, tl1);
      T cl3 = c.l + cl2;
      return algorithms::FastTwoSum(c.h, cl3);
    } else
      static_assert(sizeof(T) == 0,
                    "Fast and accurate modes are supported when using FMA");
  } else {
    // Not using FMA
    if constexpr (p == Mode::Fast) {
      // DWTimesDW1 in Joldes et al. (2017)
      two<T> c = algorithms::TwoProd(x.h, y.h);
      T tl1 = x.h * y.l;
      T tl2 = x.l * y.h;
      T cl2 = tl1 + tl2;
      T cl3 = c.l + cl2;
      return algorithms::FastTwoSum(c.h, cl3);
    } else
      static_assert(sizeof(T) == 0, "Only fast mode is supported without FMA");
  }
}

/// \brief Divides two double-word floating point numbers.
/// \details Proposed by Joldeş et al. (2017).
/// \param x The first double-word floating point number.
/// \param y The second double-word floating point number.
/// \tparam useFMA Whether to use FMA instructions.
template <bool useFMA, typename T>
inline two<T> div(const two<T> &x, T y) {
  // DWDivFP3 in Joldes et al. (2017)
  T th = x.h / y;
  two<T> pi = algorithms::TwoProd<T, useFMA>(th, y);
  T deltah = x.h - pi.h;
  T deltat = deltah - pi.l;
  T delta = deltat + x.l;
  T tl = delta / y;
  return algorithms::FastTwoSum(th, tl);
}

/// \brief Divides two double-word floating point numbers.
/// \details Proposed by Joldeş et al. (2017).
/// \param x The first double-word floating point number.
/// \param y The second double-word floating point number.
/// \tparam mode The mode (fast or accurate). When not using FMA, only the fast
/// mode is supported. Accurate mode requires double the amount of operations.
/// \tparam useFMA Whether to use FMA instructions.
template <Mode mode, bool useFMA, typename T>
inline two<T> div(const two<T> &x, const two<T> &y) {
  static_assert(mode == Mode::Fast || useFMA,
                "Only fast mode is supported without FMA");

  if constexpr (mode == Mode::Fast) {
    // DWDivDW2 in Joldes et al. (2017)
    T th = x.h / y.h;
    two<T> r = mul<Mode::Accurate, useFMA>(y, th);
    T pih = x.h - r.h;
    T deltal = x.l - r.l;
    T delta = pih + deltal;
    T tl = delta / y.h;
    return algorithms::FastTwoSum(th, tl);
  } else if constexpr (mode == Mode::Accurate) {
    // DWDivDW3 in Joldes et al. (2017)
    static_assert(useFMA, "Accurate mode requires FMA");
    T th = 1 / y.h;
    T rh = 1 - y.h * th;
    T rl = -(y.l * th);
    two<T> e = algorithms::FastTwoSum(rh, rl);
    two<T> delta = mul<Mode::Accurate, true>(e, th);
    two<T> m = add(delta, th);
    return mul<Mode::Fast, true>(x, m);
  } else
    static_assert(sizeof(T) == 0, "Unsupported mode");
}

// Reference: QD / dd_const.cpp / inline.h
// TODO: transform this to FP32 (this is FP64)
// TODO: check ranges (i.e., exponents)
const two<float> fp32_2pi  = two<float>(6.2831853e+00, 2.4492935e-16);
const two<float> fp32_pi2  = two<float>(1.5707963e+00, 6.1232339e-17);
const two<float> fp32_pi16 = two<float>(1.9634954e-01, 7.6540424e-18);
static const float fp32_nan = std::numeric_limits<float>::quiet_NaN();
const float fp32_eps = 4.9303806e-32;  // 2^-104
const double fp32_qd_split_thresh = 6.6969287e+299 // = 2^996

// TODO: check it was correctly copied from QD
// TODO: check ranges (i.e., exponents)
const two<double> fp64_2pi  = two<double>(6.283185307179586232e+00, 2.449293598294706414e-16);
const two<double> fp64_pi2  = two<double>(1.570796326794896558e+00, 6.123233995736766036e-17);
const two<double> fp64_pi16 = two<double>(1.963495408493620697e-01, 7.654042494670957545e-18);
static const double fp64_nan = std::numeric_limits<double>::quiet_NaN();
const double fp64_eps = 4.93038065763132e-32;  // 2^-104
const double fp64_qd_split_thresh = 6.69692879491417e+299 // = 2^996

// Reference: QD / inline.h
/* Computes the nearest integer to input. */
template <typename T>
inline T nint(T input) {
  T pointfive;

  if constexpr (std::is_same_v<T, float>) {
    pointfive = 0.5f;
  } else if constexpr (std::is_same_v<T, double>) {
    pointfive = 0.5;
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  if(input = std::floor(input)) {
    return input;
  }
  return std::floor(input + pointfive);
}

// Reference: QD / dd_inline.h
/* Cast to double. */
template <typename T>
inline T to_double(const two<T> &input) {
  return input.h;
}

// Reference: QD / inline.h
/* Computes high word and low word of input */
template <typename T>
inline void split(T input, T &hi, T &lo) {
  T temp;

  T local_qd_split_thresh;
  T local_split_factor;
  if constexpr (std::is_same_v<T, float>) {
    local_qd_split_thresh = fp32_qd_split_thresh;
    local_split_factor = 3.7252902e-09;  // 2^-28
  } else if constexpr (std::is_same_v<T, double>) {
    local_qd_split_thresh = fp64_qd_split_thresh;
    local_split_factor = 3.7252902984619140625e-09;  // 2^-28
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  const T qd_splitter = 134217729.0;               // = 2^27 + 1

  if(input > local_qd_split_thresh || input < -local_qd_split_thresh) {
    input *= local_split_factor; // 2^-28
    temp = qd_splitter * input;
    hi = temp - (temp - input);
    lo = input - hi;
    hi *= 268435456.0;  // 2^28
    lo *= 268435456.0;  // 2^28
  } else {
    temp = qd_splitter * input;
    hi = temp - (temp - input);
    lo = input - hi;
  }
}

// Reference: QD / inline.h
/* Computes fl(input * input) and err(input * input) */
template <typename T>
inline T two_sqr(T input, T &err) {

  T twopointzero;
  if constexpr (std::is_same_v<T, float>) {
    twopointzero = 2.0f;
  } else if constexpr (std::is_same_v<T, double>) {
    twopointzero = 2.0;
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  T hi, lo;
  T q = input * input;
  split(input, hi, lo);
  err = ((hi * hi - q) + twopointzero * hi * lo) + lo * lo;
  return q;
}

// Reference: QD / inline.h
/* Computes fl(a+b) and err(a+b).
   Assumes |a| >= |b|. */
template <typename T>
inline T quick_two_sum(T a, T b, T &err) {
  T s = a + b;
  err = b - (s - a);
  return s;
}

// Reference: QD / dd_inline.h
template <typename T>
inline two<T> sqr(const two<T> &input) {

  T twopointzero;
  if constexpr (std::is_same_v<T, float>) {
    twopointzero = 2.0f;
  } else if constexpr (std::is_same_v<T, double>) {
    twopointzero = 2.0;
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  T p1, p2;
  T s1, s2;
  p1 = two_sqr(input.h, p2);
  p2 += twopointzero * input.h * input.l;
  p2 += input.l * input.l;
  s1 = quick_two_sum(p1, p2, s2);
  return two<T>(s1, s2);
}

// Reference: QD / dd_real.cpp
/* Computes sin(a) using Taylor series.
   Assumes |a| <= pi/32 */
template<typename T>
static two<T> sin_taylor(const two<T> &input) {

  T pointfive, local_eps;
  if constexpr (std::is_same_v<T, float>) {
    pointfive = 0.5f;
    local_eps = fp32_eps;
  } else if constexpr (std::is_same_v<T, double>) {
    pointfive = 0.5;
    local_eps = fp64_eps;
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  const T thresh = pointfive * std::abs(to_double(input)) * local_eps;

  two<T> r, s, t, x;

  // TODO: make sure this is already supported in twofloat
  if (input.is_zero()) {
    return 0.0;
  }

  int i = 0;
  x = -sqr(input);

}

// Reference: QD / dd_real.cpp
template <typename T>
inline two<T> sin(const two<T> &input) {

  // TODO: make sure this is already supported in twofloat
  if (input.is_zero()) {
    return 0.0;
  }

  T local_2pi, local_pi2, pointfive, local_pi16, local_nan;
  if constexpr (std::is_same_v<T, float>) {
    local_2pi  = fp32_2pi;
    local_pi2  = fp32_pi2;
    pointfive  = 0.5f;
    local_pi16 = fp32_pi16:
    local_nan  = fp32_nan;
  } else if constexpr (std::is_same_v<T, double>) {
    local_2pi  = fp64_2pi;
    local_pi2  = fp64_pi2;
    pointfive  = 0.5;
    local_pi16 = fp64_pi16;
    local_nan  = fp64_nan;
  } else {
    std::error("LSV: other types are unsupported"); // TODO: make sure std::error is best way to proceed here
  }

  // Approximately reducing modulo 2*pi
  two<T> z = nint(input / local_2pi);
  two<T> r = input - mul(local_2pi, z);

  // Approximately reducing modulo pi/2 and then modulo pi/16

  //TODO: original type is double, here it is templated
  T q = std::floor(r.h / local_pi2.h + pointfive);
  two<T> t = r - mul(local_pi2, q);
  int j = static_cast<int>(q);
  q = std::floor(t.h / local_pi16.h + pointfive);
  t -= mul(local_pi16, q);
  int k = static_cast<int>(q);
  int abs_k = std:abs(k);

  if (j < -2 || j > 2) {
    std::error("LSV: cannot reduce modulo pi/2");
    return local_nan;
  }

  if (abs_k > 4) {
    std::error("LSV: cannot reduce modulo pi/16");
    return local_nan;
  }

}
}  // namespace doubleword
}  // namespace twofloat