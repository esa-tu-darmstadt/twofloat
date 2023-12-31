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
}  // namespace doubleword
}  // namespace twofloat