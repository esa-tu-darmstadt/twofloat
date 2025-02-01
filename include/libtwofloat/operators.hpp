#pragma once

#include <libtwofloat/arithmetics/double-word-arithmetic.hpp>

namespace twofloat {
/*********** Additions ************/
template <typename T>
inline two<T> operator+(two<T> a, two<T> b) {
  return doubleword::add<doubleword::Mode::Accurate>(a, b);
}
template <typename T>
inline two<T> operator+(T a, two<T> b) {
  return doubleword::add(a, b);
}
template <typename T>
inline two<T> operator+(two<T> a, T b) {
  return doubleword::add(a, b);
}

/*********** Subtractions ************/
template <typename T>
inline two<T> operator-(two<T> a, two<T> b) {
  return doubleword::sub<doubleword::Mode::Accurate>(a, b);
}
template <typename T>
inline two<T> operator-(T a, two<T> b) {
  return doubleword::sub(a, b);
}
template <typename T>
inline two<T> operator-(two<T> a, T b) {
  return doubleword::sub(a, b);
}

/*********** Multiplications ************/
template <typename T>
inline two<T> operator*(two<T> a, two<T> b) {
  return doubleword::mul<doubleword::Mode::Fast, false>(a, b);
}
template <typename T>
inline two<T> operator*(T a, two<T> b) {
  return doubleword::mul<doubleword::Mode::Accurate, false>(a, b);
}
template <typename T>
inline two<T> operator*(two<T> a, T b) {
  return doubleword::mul<doubleword::Mode::Accurate, false>(a, b);
}

/*********** Divisions ************/
template <typename T>
inline two<T> operator/(two<T> a, two<T> b) {
  return doubleword::div<doubleword::Mode::Fast, false>(a, b);
}
template <typename T>
inline two<T> operator/(T a, two<T> b) {
  // Maybe there is a faster implementation?
  return two<T>(a) / b;
}
template <typename T>
inline two<T> operator/(two<T> a, T b) {
  return doubleword::div<false>(a, b);
}

}  // namespace twofloat