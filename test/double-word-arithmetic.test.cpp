#include <libtwofloat/arithmetics/double-word-arithmetic.hpp>
#include <libtwofloat/limits.hpp>
#include <limits>
#include <random>

#include "gtest/gtest.h"

using namespace twofloat;

typedef two<float> floatfloat;

namespace twofloat {
namespace doubleword {
namespace test {

template <doubleword::Mode mode>
void addDWTest() {
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
  two<float> result = doubleword::add<mode>(x, y);

  // Calculate the result using single precision
  float fresult = (float)dx + (float)dy;

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

template <doubleword::Mode mode>
void subDWTest() {
  // This test subtracts two numbers with different exponents, which normally
  // rounds away information stored in the lower bits of the smaller number's
  // mantissa.
  double dx = 10;
  double dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx - dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::sub<mode>(x, y);

  // Calculate the result using single precision
  float fresult = (float)dx - (float)dy;

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

template <doubleword::Mode mode, bool useFMA>
void mulDWTest() {
  // This test multiplies two numbers with information stored in the lower bits
  // of the mantissa. During a normal float multiplication, this information is
  // lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  float dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx * dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::mul<mode, useFMA>(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

template <doubleword::Mode mode, bool useFMA>
void mulFPTest() {
  // This test multiplies two numbers with information stored in the lower bits
  // of the mantissa. During a normal float multiplication, this information is
  // lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  float dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx * dy;

  two<float> x = {(float)dx, 0.0f};
  float y = (float)dy;

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::mul<mode, useFMA>(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

template <doubleword::Mode mode, bool useFMA>
void divDWTest() {
  // These are two numbers that set the last bit of the mantissa to 1.
  // During a normal float multiplication, this information is lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  float dy = 10;

  // Calculate the result using double precision
  double dresult = dx / dy;

  two<float> x = {(float)dx, 0.0f};
  two<float> y = {(float)dy, 0.0f};

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::div<mode, useFMA>(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

template <bool useFMA>
void divFPTest() {
  // These are two numbers that set the last bit of the mantissa to 1.
  // During a normal float multiplication, this information is lost.
  double dx = 1.0 + std::numeric_limits<float>::epsilon();
  float dy = 10;

  // Calculate the result using double precision
  double dresult = dx / dy;

  two<float> x = {(float)dx, 0.0f};
  float y = (float)dy;

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::div<useFMA>(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

TEST(DoubleWordArithmetic, AddDWAccurateTest) {
  addDWTest<doubleword::Mode::Accurate>();
}

TEST(DoubleWordArithmetic, AddDWSloppyTest) {
  addDWTest<doubleword::Mode::Sloppy>();
}

TEST(DoubleWordArithmetic, AddFPTest) {
  // This test adds two numbers with different exponents, which normally
  // rounds away information stored in the lower bits of the smaller number's
  // mantissa.
  double dx = 10;
  double dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx + dy;

  two<float> x = {(float)dx, 0.0f};
  float y = (float)dy;

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::add(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

TEST(DoubleWordArithmetic, SubDWAccurateTest) {
  subDWTest<doubleword::Mode::Accurate>();
}

TEST(DoubleWordArithmetic, SubDWSloppyTest) {
  subDWTest<doubleword::Mode::Sloppy>();
}

TEST(DoubleWordArithmetic, SubFPTest) {
  // This test adds two numbers with different exponents, which normally
  // rounds away information stored in the lower bits of the smaller number's
  // mantissa.
  double dx = 10;
  double dy = 1.0 + std::numeric_limits<float>::epsilon();

  // Calculate the result using double precision
  double dresult = dx - dy;

  two<float> x = {(float)dx, 0.0f};
  float y = (float)dy;

  // Calculate the result using pair arithmetic
  two<float> result = doubleword::sub(x, y);

  // Expect the double-word result to be close to the double result,
  // considering the precision of the double-word arithmetic.
  EXPECT_NEAR(dresult, result.eval<double>(),
              dresult * std::numeric_limits<two<float>>::epsilon());
}

TEST(DoubleWordArithmetic, MulDWAccurateTest) {
  mulDWTest<doubleword::Mode::Fast, false>();
}

TEST(DoubleWordArithmetic, MulDWFMAFastTest) {
  mulDWTest<doubleword::Mode::Fast, true>();
}

TEST(DoubleWordArithmetic, MulDWFMAAccurateTest) {
  mulDWTest<doubleword::Mode::Accurate, true>();
}

TEST(DoubleWordArithmetic, MulFPAccurateTest) {
  mulFPTest<doubleword::Mode::Fast, false>();
}

TEST(DoubleWordArithmetic, MulFPFMAFastTest) {
  mulFPTest<doubleword::Mode::Fast, true>();
}

TEST(DoubleWordArithmetic, MulFPFMAAccurateTest) {
  mulFPTest<doubleword::Mode::Accurate, true>();
}

TEST(DoubleWordArithmetic, DivDWAccurateTest) {
  divDWTest<doubleword::Mode::Fast, false>();
}

TEST(DoubleWordArithmetic, DivDWFMAFastTest) {
  divDWTest<doubleword::Mode::Fast, true>();
}

TEST(DoubleWordArithmetic, DivDWFMAAccurateTest) {
  divDWTest<doubleword::Mode::Accurate, true>();
}

TEST(DoubleWordArithmetic, DivFPTest) { divFPTest<false>(); }

TEST(DoubleWordArithmetic, DivFPFMATest) { divFPTest<true>(); }

TEST(DoubleWordArithmetic, SinTest) {
  float x = 1;
  float sinx = std::sin(x);

  two<float> two_x = {x, 0.0f};
  two<float> sin_two_x = doubleword::sin(two_x);

  // Check the expected and actual values within /build/Testing/Temporary/LastTest.log
  EXPECT_NEAR(sinx, sin_two_x.eval<float>(), sinx * std::numeric_limits<two<float>>::epsilon());
}

TEST(DoubleWordArithmetic, SinTest2) {
  double x = 1;
  double sinx = std::sin(x);

  two<double> two_x = {x, 0.0};
  two<double> sin_two_x = doubleword::sin(two_x);

  // Check the expected and actual values within /build/Testing/Temporary/LastTest.log
  EXPECT_NEAR(sinx, sin_two_x.eval<double>(), sinx * std::numeric_limits<two<double>>::epsilon());
}

}  // namespace test
}  // namespace doubleword
}  // namespace twofloat