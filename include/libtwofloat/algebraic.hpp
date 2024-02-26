#pragma once

#include <cmath>
#include <libtwofloat/algorithms.hpp>
#include <libtwofloat/twofloat.hpp>

namespace twofloat {

template <typename T>
two<T> nint(two<T> a) {
  T hi = round(a.h);
  if (hi == a.h) {
    // high word is an integer already. Round the low word
    T lo = round(a.l);

    // renormalize and return
    return algorithms::FastTwoSum(hi, lo);
  } else {
    T lo = 0.0;
    if (abs(hi - a.h) == 0.5 && a.l < 0.0) {
      // there is a tie in the high word, consult the low word to break the tie
      hi -= 1.0;
    }
    return {hi, lo};
  }
}
}  // namespace twofloat