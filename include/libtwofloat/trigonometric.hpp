#pragma once

#include <math.h>

#include <array>
#include <cmath>
#include <libtwofloat/algebraic.hpp>
#include <libtwofloat/algorithms.hpp>
#include <libtwofloat/twofloat.hpp>

namespace twofloat {

template <typename T>
using ArithFunctionDWDW = two<T> (*)(const two<T>&, const two<T>&);
template <typename T>
using ArithFunctionDWF = two<T> (*)(const two<T>&, T);

namespace details {
template <typename T>
two<T> m_2_pi = algorithms::Split<T, long double>(
    6.283185307179586476925286766559005768394338798750211641949L);
template <typename T>
two<T> m_pi = algorithms::Split<T, long double>(
    3.1415926535897932384626433832795028841971693993751058209745L);
template <typename T>
two<T> m_pi_2 = algorithms::Split<T, long double>(
    1.5707963267948966192313216916397514420985846996875529104874L);
template <typename T>
two<T> m_pi_16 = algorithms::Split<T, long double>(
    0.1963495408493620774039152114549689302623230874609441138109L);

template <typename T>
const std::array<two<T>, 12> SinCoefficients = {
    algorithms::Split<T, long double>(
        0.16666666666666666666666666666666666606732416116558L),
    algorithms::Split<T, long double>(
        0.0083333333333333333333333333333331135404851288270047L),
    algorithms::Split<T, long double>(
        0.00019841269841269841269841269839935785325638310428717L),
    algorithms::Split<T, long double>(
        0.27557319223985890652557316053039946268333231205686e-5L),
    algorithms::Split<T, long double>(
        0.25052108385441718775048214826384312253862930064745e-7L),
    algorithms::Split<T, long double>(
        0.16059043836821614596571832194524392581082444805729e-9L),
    algorithms::Split<T, long double>(
        0.76471637318198151807063387954939213287488216303768e-12L),
    algorithms::Split<T, long double>(
        0.28114572543451292625024967174638477283187397621303e-14L),
    algorithms::Split<T, double>(
        0.82206352458348947812512122163446202498005154296863e-17L),
    algorithms::Split<T, double>(
        0.19572940011906109418080609928334380560135358385256e-19L),
    algorithms::Split<T, double>(
        0.38680813379701966970673724299207480965452616911420e-22L),
    algorithms::Split<T, double>(
        0.64038150078671872796678569586315881020659912139412e-25L)};

template <typename T, ArithFunctionDWDW<T> mul, ArithFunctionDWDW<T> add>
two<T> sin_taylor(const two<T>& input) {
  if (input.eval() == 0) return {0, 0};

  two<T> x = -(mul(input, input));
  two<T> s = input;
  two<T> p = input;

  for (int i = 0; i < SinCoefficients<T>.size(); i++) {
    p = mul(p, x);
    two<T> t = mul(p, SinCoefficients<T>[i]);
    s = add(s, t);
  }

  return s;
}
}  // namespace details

template <typename T, ArithFunctionDWDW<T> mulDWDW, ArithFunctionDWF<T> mulDWF,
          ArithFunctionDWDW<T> addDWDW, ArithFunctionDWDW<T> subDWDW,
          ArithFunctionDWDW<T> divDWDW>
two<T> sin(const two<T>& a) {
  if (a.eval() == 0) return {0, 0};

  // approximately reduce modulo 2*pi
  two<T> z = nint(divDWDW(a, details::m_2_pi<T>));
  two<T> r = subDWDW(a, mulDWDW(details::m_2_pi<T>, z));

  // approximately reduce modulo pi/2 and then modulo pi/16.
  // t = a % (pi/16)
  T q = floor(r.h / details::m_pi_2<T>.h + 0.5);
  two<T> t = subDWDW(r, mulDWF(details::m_pi_2<T>, q));
  int j = (int)q;
  q = floor(r.h / details::m_pi_16<T>.h + 0.5);
  t = subDWDW(t, mulDWF(details::m_pi_16<T>, q));
  int k = (int)q;
  int abs_k = abs(k);

  if (j < -2 || j > 2) {
    // cannot reduce modulo pi/2
    return {0, 0};
  }

  if (abs_k > 4) {
    // cannot reduce module pi/16
    return {0, 0};
  }

  // FIXME: Implement cos_taylor
  if (k == 0) {
    switch (j) {
      case 0:
        return details::sin_taylor<T, mulDWDW, addDWDW>(t);
      // case 1:
      //   return details::cos_taylor<t, mul, add>(t);
      // case -1:
      //   return details::cos_taylor(t);
      default:
        return -details::sin_taylor<T, mulDWDW, addDWDW>(t);
    }
  }

  // FIXME: The following code requires the SinTable and CosTable coefficients
  // and the sincos_taylor function to be implemented.

  // two<T> u = details::CosTable<T>[abs_k - 1];
  // two<T> v = details::SinTable<T>[abs_k - 1];

  // auto [sin_t, cos_t] = details::sincos_taylor(t);
  // if (j == 0) {
  //   if (k > 0)
  //     return add(mulDWDW(u, sin_t), mulDWDW(v, cos_t));
  //   else
  //     return sub(mulDWDW(u, sin_t), mulDWDW(v, cos_t));
  // } else if (j == 1) {
  //   if (k > 0)
  //     return sub(mulDWDW(u, cos_t), mulDWDW(v, sin_t));
  //   else
  //     return add(mulDWDW(u, cos_t), mulDWDW(v, sin_t));
  // } else {
  //   if (k > 0)
  //     return sub(-mulDWDW(u, sin_t), mulDWDW(v, cos_t));
  //   else
  //     return sub(mulDWDW(v, cos_t), mulDWDW(u, sin_t));
  // }
}

}  // namespace twofloat