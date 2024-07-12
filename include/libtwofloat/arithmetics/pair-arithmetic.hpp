#pragma once

/// \file pair-arithmetic.hpp
/// \brief Implements the pairwise arithmetic proposed by Lange and Rump.

#include <libtwofloat/algorithms.hpp>
#include <libtwofloat/twofloat.hpp>

namespace twofloat {
/// \brief Implements the pairwise arithmetic proposed by Lange and Rump in
/// 2020.
/// \details Pairwise arithmetic represents a floating point number as a
/// sum of two floating point numbers. In contrast to double-word (or
/// double-double) arithmetic, pairwise arithmetic skips the normalization step
/// after each operation. The downside of this is that the arithmetic only
/// garantuees faithful rounding up to a certain number of operations (e.g., the
/// number of summations or multiplications executed in any order). The bound on
/// the number of operations is presented in Collary 6.1 - 6.5 for
/// multiplications, polynomial evaluations, summations, dot products and the
/// the euclidian norm respectively. See https://doi.org/10.1145/3290955 for
/// details.
namespace pair {

/// \brief Adds two pairwise floating point numbers using the pairwise
/// arithmetic. This is algorithm `CPairSum` in chapter 3 of the paper.
template <typename T>
inline two<T> add(const two<T> &x, const two<T> &y) {
  two<T> s = algorithms::TwoSum(x.h, y.h);
  float v = x.l + y.l;
  float w = s.l + v;
  return {s.h, w};
}

/// \brief Adds a pairwise floating point number with a normal floating point
/// number using the pairwise arithmetic. This is derived from algorithm
/// `CPairSum` in chapter 3 of the paper.
template <typename T>
inline two<T> add(const two<T> &x, T y) {
  two<T> s = algorithms::TwoSum(x.h, y);
  float w = s.l + x.l;
  return {s.h, w};
}
template <typename T>
inline two<T> add(T x, const two<T> &y) {
  return add(y, x);
}

/// \brief Subtracts two pairwise floating point numbers using the pairwise
/// arithmetic. This is derived from algorithm `CPairDiff` in chapter 3 of the
/// paper.
template <typename T>
inline two<T> sub(const two<T> &x, const two<T> &y) {
  two<T> s = algorithms::TwoDiff(x.h, y.h);
  float v = x.l - y.l;
  float w = s.l + v;
  return {s.h, w};
}

/// \brief Subtracts a pairwise floating point number with a normal floating
/// point number using the pairwise arithmetic. This is derived from algorithm
/// `CPairDiff` in chapter 3 of the paper.
template <typename T>
inline two<T> sub(const two<T> &x, T y) {
  two<T> s = algorithms::TwoDiff(x.h, y);
  float w = s.l + x.l;
  return {s.h, w};
}
/// \brief Subtracts a pairwise floating point number with a normal floating
/// point number using the pairwise arithmetic. This is derived from algorithm
/// `CPairDiff` in chapter 3 of the paper.
template <typename T>
inline two<T> sub(T x, const two<T> &y) {
  two<T> s = algorithms::TwoDiff(x, y.h);
  float w = s.l + y.l;
  return {s.h, w};
}

/// \brief Multiplies two pairwise floating point numbers using the pairwise
/// arithmetic. This is algorithm `CPairProd` in chapter 3 of the paper.
template <bool useFMA, typename T>
inline two<T> mul(const two<T> &x, const two<T> &y) {
  two<T> c = algorithms::TwoProd<T, useFMA>(x.h, y.h);
  T tl1 = x.h * y.l;
  T tl2 = x.l * y.h;
  T cl2 = tl1 + tl2;
  T cl3 = c.l + cl2;
  return {c.h, cl3};
}

/// \brief Multiplies a pairwise floating point number with a normal floating
/// point number using the pairwise arithmetic. This is derived from algorithm
/// `CPairProd` in chapter 3 of the paper.
template <bool useFMA, typename T>
inline two<T> mul(const two<T> &x, T y) {
  two<T> c = algorithms::TwoProd<T, useFMA>(x.h, y);
  T tl2 = x.l * y;
  T cl3 = c.l + tl2;
  return {c.h, cl3};
}
template <bool useFMA, typename T>
inline two<T> mul(T x, const two<T> &y) {
  return mul<useFMA>(y, x);
}

/// \brief Multiplies two normal floating point number using the pairwise
/// arithmetic. This is derived from algorithm `CPairProd` in chapter 3 of the
/// paper.
template <bool useFMA, typename T>
inline two<T> mul(T x, T y) {
  return algorithms::TwoProd<T, useFMA>(x, y);
}

/// \brief Divides two pairwise floating point numbers using the pairwise
/// arithmetic. This is algorithm `CPairDiv` in chapter 3 of the paper.
template <typename T>
inline two<T> div(const two<T> &x, const two<T> &y) {
  T c = x.h / y.h;
  T t = x.h - c * y.h;
  T p = t + x.l;
  T q = c * y.l;
  T r = p - q;
  T s = y.h + y.l;
  T g = r / s;
  return {c, g};
}

/// \brief Divides a pairwise floating point number by a floating point number.
/// This is derived from algorithm `CPairDiv` in chapter 3 of the
/// paper.
template <typename T>
inline two<T> div(const two<T> &x, T y) {
  T c = x.h / y;
  T t = x.h - c * y;
  T p = t + x.l;
  T r = p;
  T s = y;
  T g = r / s;
  return {c, g};
}

}  // namespace pair
}  // namespace twofloat
