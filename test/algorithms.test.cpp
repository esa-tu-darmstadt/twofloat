#include <libtwofloat/algorithms.hpp>

#include <cmath>

#include "gtest/gtest.h"

namespace twofloat {
namespace algorithms {
namespace test {

TEST(Algorithms, FmaPreservesTemplatePrecision) {
  double a = std::nextafter(1.0, 2.0);

  EXPECT_EQ(std::fma(a, a, -1.0), algorithms::fma(a, a, -1.0));
}

}  // namespace test
}  // namespace algorithms
}  // namespace twofloat
