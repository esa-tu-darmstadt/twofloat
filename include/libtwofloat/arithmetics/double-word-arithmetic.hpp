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

/// \brief QD constants
/// Reference: QD / dd_const.cpp / inline.h / dd_real.cpp
template <typename T>
static const two<T> _2pi = {6.283185307179586232e+00, 2.449293598294706414e-16};

template <typename T>
static const two<T> _pi2  = {1.570796326794896558e+00, 6.123233995736766036e-17};

template <typename T>
static const two<T> _pi16 = {1.963495408493620697e-01, 7.654042494670957545e-18};

template <typename T>
static const constexpr T _nan = std::numeric_limits<T>::quiet_NaN();

template <typename T>
static const constexpr T _eps = 4.93038065763132e-32;  // 2^-104

static const int n_inv_fact = 15;

/// \brief Some recurrently-used contants
template <typename T>
static const constexpr T pointfive = 0.5;

template <typename T>
static const constexpr T twopointzero = 2.0;

template <typename T>
static const constexpr T zeropointzero = 0.0;

template <typename T>
static const constexpr T onepointzero = 1.0;

/// \brief Tables containing constant data for trigonometric functions
template <typename T>
struct constants_trig_tables {

  static const constexpr T inv_fact[n_inv_fact][2] = {
    { 1.66666666666666657e-01,  9.25185853854297066e-18},
    { 4.16666666666666644e-02,  2.31296463463574266e-18},
    { 8.33333333333333322e-03,  1.15648231731787138e-19},
    { 1.38888888888888894e-03, -5.30054395437357706e-20},
    { 1.98412698412698413e-04,  1.72095582934207053e-22},
    { 2.48015873015873016e-05,  2.15119478667758816e-23},
    { 2.75573192239858925e-06, -1.85839327404647208e-22},
    { 2.75573192239858883e-07,  2.37677146222502973e-23},
    { 2.50521083854417202e-08, -1.44881407093591197e-24},
    { 2.08767569878681002e-09, -1.20734505911325997e-25},
    { 1.60590438368216133e-10,  1.25852945887520981e-26},
    { 1.14707455977297245e-11,  2.06555127528307454e-28},
    { 7.64716373181981641e-13,  7.03872877733453001e-30},
    { 4.77947733238738525e-14,  4.39920548583408126e-31},
    { 2.81145725434552060e-15,  1.65088427308614326e-31}
  };

  // Reference: QD / dd_real.cpp
  // Table of sin(k * pi/16) and cos(k * pi/16).
  static const constexpr T cos_table [4][2] = {
    {9.807852804032304306e-01, 1.854693999782500573e-17},
    {9.238795325112867385e-01, 1.764504708433667706e-17},
    {8.314696123025452357e-01, 1.407385698472802389e-18},
    {7.071067811865475727e-01, -4.833646656726456726e-17}
  };
  static const constexpr T sin_table [4][2] = {
    {1.950903220161282758e-01, -7.991079068461731263e-18},
    {3.826834323650897818e-01, -1.005077269646158761e-17},
    {5.555702330196021776e-01,  4.709410940561676821e-17},
    {7.071067811865475727e-01, -4.833646656726456726e-17}
  };
};

/// \brief Computes fl(input * input) and err(input * input)
/// Reference: QD / inline.h (within qd namespace)
template <typename T>
inline T two_sqr(T input, T &err) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");
  T hi, lo;
  T q = input * input;
  two<T> temp = twofloat::algorithms::Split(input);
  hi = temp.h;
  lo = temp.l;
  err = ((hi * hi - q) + twopointzero<T> * hi * lo) + lo * lo;
  return q;
}

/// Reference QD / inline.h (within qd namespace)
/// TODO: make sure this has not been already implemented in twofloat
template <typename T>
inline T two_sum(T a, T b, T &err) {
  T s = a + b;
  T bb = s - a;
  err = (a - (s - bb)) - (b - bb);
  return s;
}

/// Reference: QD / dd_inline.h (within dd_real namespace)
/// TODO: make sure this has not been already implemented in twofloat
template <typename T>
inline two<T> add(T a, T b) {
  T s, e;
  s = two_sum(a, b, e);
  return two<T> (s, e);
}

/// \brief Computes the nearest integer to input.
/// Reference: QD / inline.h (within qd namespace)
template <typename T>
inline T nint(T input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");
  if(input = std::floor(input)) {
    return input;
  }
  return std::floor(input + pointfive<T>);
}

/// \brief Round to nearest integer.
/// Reference: QD / dd_inline.h
template <typename T>
inline two<T> nint(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");
  T hi = nint(input.h);
  T lo;
  if(hi == input.h) {
    // High word is an integer already. Round the low word
    lo = nint(input.l);

    // Renormalize. This is needed if h = some integer, l = 1/2
    //hi = qd::quick_two_sum(hi, lo, lo);
    two<T> temp = twofloat::algorithms::FastTwoSum(hi, lo);
    hi = temp.h;
    lo = temp.l;
  }
  else {
    // High word is not an integer
    lo = zeropointzero<T>;
    if(std::abs(hi - input.h) == pointfive<T> && input.l < zeropointzero<T>) {
      // There is a tie in the high word, consult the low word to break the tie
      hi -= onepointzero<T>; /* NOTE: This does not cause INEXACT. */
    }
  }

  return two<T>{hi, lo};
}

/// Reference: QD / dd_inline.h (within dd_real namespace)
template <typename T>
inline two<T> sqr(T input) {
  T p1, p2;
  p1 = two_sqr(input, p2);
  return two<T>{p1, p2};
}

/// Reference: QD / dd_inline.h
template <typename T>
inline two<T> sqr(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");
  T p1, p2;
  T s1, s2;
  p1 = two_sqr(input.h, p2);
  p2 += twopointzero<T> * input.h * input.l;
  p2 += input.l * input.l;
  two<T> temp = twofloat::algorithms::FastTwoSum(p1, p2);
  s1 = temp.h;
  s2 = temp.l;
  return two<T>{s1, s2};
}

/// \brief Computes sin(a) using Taylor series.
/// Assumes |a| <= pi/32
/// Reference: QD / dd_real.cpp
template <Mode p, bool useFMA, typename T>
static two<T> sin_taylor(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");

  const T* local_ptr_inv_fact = &constants_trig_tables<T>::inv_fact[0][0];
  T local_eps = _eps<T>;
  const T thresh = pointfive<T> * std::abs(input.eval()) * local_eps;
  two<T> r, s, t, x;

  if (input.eval() == zeropointzero<T>) {
    return two<T>{zeropointzero<T>, zeropointzero<T>};
  }

  int i = 0;
  two<T> temp = sqr(input);
  x = two<T>{-temp.h, -temp.l};
  s = input;
  r = input;

  do {
    r = mul<p, useFMA>(r, x);
    t = mul<p, useFMA>(r, two<T>{(local_ptr_inv_fact+i)[0], (local_ptr_inv_fact+i)[1]});
    s = add<p>(s, t);
    i += 2;
  } while (i < n_inv_fact && std::abs(t.eval()) > thresh);

  return s;
}

/// \brief double-double * double, where double is a power of 2.
/// Reference: QD / dd_inline.h
template<typename T>
inline two<T> mul_pwr2(const two<T> &input, T b) {
  return two<T>(input.h * b, input.l * b);
}

/// Reference: QD / dd_real.cpp
template <Mode p, bool useFMA, typename T>
static two<T> cos_taylor(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");

  const T* local_ptr_inv_fact = &constants_trig_tables<T>::inv_fact[0][0];
  T local_eps = _eps<T>;
  const T thresh = pointfive<T> * local_eps;
  two<T> r, s, t, x;

  if (input.eval() == zeropointzero<T>) {
    return two<T>{onepointzero<T>, zeropointzero<T>};
  }

  two<T> temp = sqr(input);
  x = two<T>{-temp.h, -temp.l};
  r = x;
  s = add(onepointzero<T>, mul_pwr2(r, pointfive<T>));
  int i = 1;
  do {
    r = mul<p, useFMA>(r, x);
    t = mul<p, useFMA>(r, two<T>{(local_ptr_inv_fact+i)[0], (local_ptr_inv_fact+i)[1]});
    s = add<p>(s, t);
    i += 2;
  } while (i < n_inv_fact && std::abs(t.eval()) > thresh);

  return s;
}

/// \brief Computes the square root of the double-double number dd.
/// NOTE: dd must be a non-negative number.
///
/// Strategy:  Use Karp's trick:  if x is an approximation
///     to sqrt(a), then
///
///        sqrt(a) = a*x + [a - (a*x)^2] * x / 2   (approx)
///
///     The approximation is accurate to twice the accuracy of x.
///     Also, the multiplication (a*x) and [-]*x can be done with
///     only half the precision.
///
/// Reference: QD / dd_real.cpp
template <Mode p, typename T>
two<T> sqrt(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");
  T local_nan = _nan<T>;

  if (input.eval() == zeropointzero<T>) {
    return two<T>{zeropointzero<T>, zeropointzero<T>};
  }

  // ERROR: Negative argument
  // if (input.is_negative())
  if (input.h < zeropointzero<T>) {
    return two<T>{local_nan, local_nan};
  }

  T x = onepointzero<T> / std::sqrt(input.h);
  T ax = input.h * x;
  two<T> temp = sub<p>(input, sqr(ax));
  return add(ax, temp.h * x * pointfive<T>);
}

/// Reference: QD / dd_real.cpp
template <Mode p, bool useFMA, typename T>
static void sincos_taylor(const two<T> &input, two<T> &sin_a, two<T> &cos_a) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");

  if (input.eval() == zeropointzero<T>) {
    sin_a.h = zeropointzero<T>;
    sin_a.l = zeropointzero<T>;
    cos_a.h = onepointzero<T>;
    cos_a.l = zeropointzero<T>;
    return;
  }

  sin_a = sin_taylor<p, useFMA>(input);
  cos_a = sqrt<p>(sub(onepointzero<T>, sqr(sin_a)));
}

/// Reference: QD / dd_real.cpp
template <Mode p, bool useFMA, typename T>
inline two<T> sin(const two<T> &input) {
  static_assert((std::is_same_v<T, float> || std::is_same_v<T, double>), "Other types not supported");

  two<T> local_2pi = _2pi<T>;
  two<T> local_pi2 = _pi2<T>;
  two<T> local_pi16 = _pi16<T>;
  T local_nan = _nan<T>;

  const T* local_ptr_cos_table = &constants_trig_tables<T>::cos_table[0][0];
  const T* local_ptr_sin_table = &constants_trig_tables<T>::sin_table[0][0];

  if (input.eval() == zeropointzero<T>) {
    return two<T>{zeropointzero<T>, zeropointzero<T>};
  }

  // Approximately reducing modulo 2*pi
  two<T> z = nint(div<p, useFMA>(input, local_2pi));
  two<T> r = sub<p>(input, mul<p, useFMA>(local_2pi, z));

  /*
  std::cout << "LSV starts ..." << std::endl;
  std::cout << "input = " << input.eval() << std::endl;
  std::cout << "z = " << z.eval() << std::endl;
  std::cout << "r = " << r.eval() << std::endl;
  std::cout << "... LSV ends" << std::endl;
  */

  // Approximately reducing modulo pi/2 and then modulo pi/16
  T q = std::floor(r.h / local_pi2.h + pointfive<T>);
  two<T> t = sub<p>(r, mul<p, useFMA>(local_pi2, q));
  int j = static_cast<int>(q);
  q = std::floor(t.h / local_pi16.h + pointfive<T>);
  t = sub<p>(t, mul<p, useFMA>(local_pi16, q));
  int k = static_cast<int>(q);
  int abs_k = std::abs(k);

  /*
  std::cout << "LSV starts ..." << std::endl;
  std::cout << "t = " << t.eval() << std::endl;
  std::cout << "q = " << q << std::endl;
  std::cout << "j = " << j << std::endl;
  std::cout << "k = " << k << std::endl;
  std::cout << "abs_k = " << abs_k << std::endl;
  std::cout << "... LSV ends" << std::endl;
  */

  // ERROR: cannot reduce modulo pi/2\n
  if (j < -2 || j > 2) {
    return two<T>{local_nan, local_nan};
  }

  // ERROR: cannot reduce modulo pi/16
  if (abs_k > 4) {
    return two<T>{local_nan, local_nan};
  }

  if (k == 0) {
    switch(j) {
      case 0:
        return sin_taylor<p, useFMA>(t);
      case 1:
        return cos_taylor<p, useFMA>(t);
      case -1:
        return two<T>{-cos_taylor<p, useFMA>(t).h, -cos_taylor<p, useFMA>(t).l}; // Reference: QD / dd_inline.h
      default:
        return two<T>{-sin_taylor<p, useFMA>(t).h, -sin_taylor<p, useFMA>(t).l};
    }
  }

  // IMPORTANT: notice the factor 2 when accessing the cos() and sin() look-up tables
  // This factor is required to access data correctly as done in QD
  two<T> u{(local_ptr_cos_table + 2 * (abs_k - 1))[0], (local_ptr_cos_table + 2 * (abs_k - 1))[1]};
  two<T> v{(local_ptr_sin_table + 2 * (abs_k - 1))[0], (local_ptr_sin_table + 2 * (abs_k - 1))[1]};

  two<T> sin_t, cos_t;

  sincos_taylor<p, useFMA>(t, sin_t, cos_t);

  /*
  std::cout << "LSV starts ..." << std::endl;
  std::cout << "sin_t = " << sin_t.eval() << std::endl;
  std::cout << "cos_t = " << cos_t.eval() << std::endl;
  std::cout << "... LSV ends" << std::endl;
  */

  if (j == 0) {
    if (k > 0) {
      r = add<p>(mul<p, useFMA>(u, sin_t), mul<p, useFMA>(v, cos_t));
    } else {
      r = sub<p>(mul<p, useFMA>(u, sin_t), mul<p, useFMA>(v, cos_t));
    }
  } else if (j == 1) {
    if (k > 0) {
      r = sub<p>(mul<p, useFMA>(u, cos_t), mul<p, useFMA>(v, sin_t));
    } else {
      r = add<p>(mul<p, useFMA>(u, cos_t), mul<p, useFMA>(v, sin_t));
      /*
      std::cout << "LSV starts ..." << std::endl;
      std::cout << "u = " << u.eval() << std::endl;
      std::cout << "cos_t = " << cos_t.eval() << std::endl;
      std::cout << "v = " << v.eval() << std::endl;
      std::cout << "sin_t = " << sin_t.eval() << std::endl;
      std::cout << "r = " << r.eval() << std::endl;
      std::cout << "... LSV ends" << std::endl;
      */
    }
  } else if (j == -1) {
    if (k > 0) {
      r = sub<p>(mul<p, useFMA>(v, sin_t), mul<p, useFMA>(u, cos_t));
    } else if (k < 0) {
      two<T> temp = mul<p, useFMA>(u, cos_t);
      two<T> neg_temp = two<T>{-temp.h, -temp.l};
      r = sub<p>(neg_temp, mul<p, useFMA>(v, sin_t));
    }
  } else {
    if (k > 0) {
      two<T> temp = mul<p, useFMA>(u, sin_t);
      two<T> neg_temp = two<T>{-temp.h, -temp.l};
      r = sub<p>(neg_temp, mul<p, useFMA>(v, cos_t));
    } else {
      r = sub<p>(mul<p, useFMA>(v, cos_t), mul<p, useFMA>(u, sin_t));
    }
  }
  /*
  std::cout << "LSV starts ..." << std::endl;
  std::cout << "r = " << r.eval() << std::endl;
  std::cout << "... LSV ends" << std::endl;
  */
  return r;
}
}  // namespace doubleword
}  // namespace twofloat