#include <libtwofloat/arithmetics/pair-arithmetic.hpp>
#include <libtwofloat/limits.hpp>
#include <random>

#include "gtest/gtest.h"

using namespace twofloat;

typedef two<float> floatfloat;

namespace twofloat {
namespace pair {
namespace test {

TEST(PairArithmetic, AddTest) {
  // This test adds two numbers with different exponents, which normally
  // rounds away information stored in the lower bits of the smaller number's
  // mantissa.
  double dx = 10;
  double dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx + dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = pair::add(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

TEST(PairArithmetic, SubTest) {
  // This test subtracts two numbers with different exponents, which normally
  // rounds away information stored in the lower bits of the smaller number's
  // mantissa.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  double dy = 10;

  // Calculate the result using double precision
  double dresult = dx - dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = pair::sub(x, y);

  // Calculate the result using single precision
  float fresult = (float)dx - (float)dy;

  // Expect the float results to be different from the double result
  EXPECT_NE(dresult, fresult);

  // Expect the pair results to be equal to the double result
  EXPECT_EQ(dresult, result.eval<double>());
}

TEST(PairArithmetic, MulTest) {
  // These are two numbers that set the last bit of the mantissa to 1.
  // During a normal float multiplication, this information is lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  double dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx * dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = pair::mul<false>(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

TEST(PairArithmetic, DivTest) {
  // These are two numbers that set the last bit of the mantissa to 1.
  // During a normal float multiplication, this information is lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  float dy = 10;

  // Calculate the result using double precision
  double dresult = dx / dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = pair::div(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

}  // namespace test
}  // namespace pair
}  // namespace twofloat