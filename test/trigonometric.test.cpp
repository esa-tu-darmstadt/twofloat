#include <ios>
#include <libtwofloat/arithmetics/double-word-arithmetic.hpp>
#include <libtwofloat/limits.hpp>
#include <libtwofloat/trigonometric.hpp>
#include <limits>
#include <random>

#include "gtest/gtest.h"

using namespace twofloat;

typedef two<float> floatfloat;

namespace twofloat {
namespace doubleword {
namespace test {

template <doubleword::Mode mode, bool useFMA>
void sinTest(double dx) {
  // Print SinCoefficients
  // TODO: Remove this debug ouput.
  for (int i = 0; i < details::SinCoefficients<float>.size(); i++) {
    std::cout << std::scientific << std::setprecision(17);
    std::cout << "SinCoefficients[" << i
              << "] = " << details::SinCoefficients<float>[i].h << "+"
              << details::SinCoefficients<float>[i].l << std::endl;
  }

  two<float> x = algorithms::Split<float, double>(dx);
  two<float> res = sin<float, doubleword::mul<mode, useFMA>,
                       doubleword::mul<mode, useFMA>, doubleword::add<mode>,
                       doubleword::sub<mode>, doubleword::div<mode, useFMA>>(x);

  double dresult = std::sin(dx);

  EXPECT_NEAR(dresult, res.eval<double>(), 1e-9);
}

TEST(TrigonometricTest, SinAccurateFMA) {
  for (double dx = -0.05; dx <= 0.05; dx += 0.01) {
    sinTest<doubleword::Mode::Accurate, true>(dx);
  }
}

}  // namespace test
}  // namespace doubleword
}  // namespace twofloat