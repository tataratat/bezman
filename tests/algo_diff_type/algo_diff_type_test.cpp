#include <gtest/gtest.h>
#include <cmath>

#include "bezierManipulation/src/utils/computational_derivation/algo_diff_type.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::algo_diff_type_test {

using ADT = utils::computational_derivation::AlgoDiffType<double>;

// Constructors and Basic Operations
TEST(AlgoTypeTest, TestValueCorrectnessBasic) {
  // Define some Variables
  ADT x{3., 2, 0};  // x = 3
  ADT y{2., 2, 1};  // y = 2

  // Test Values
  EXPECT_FLOAT_EQ((x + y).GetValue(), 3. + 2.);
  EXPECT_FLOAT_EQ((x * y).GetValue(), 3. * 2.);
  EXPECT_FLOAT_EQ((x / y).GetValue(), 3. / 2.);
  EXPECT_FLOAT_EQ((x - y).GetValue(), 3. - 2.);
};

TEST(AlgoTypeTest, TestDerivCorrectnessBasic) {
  // Define some Variables
  const ADT x{3., 2, 0};  // x = 3
  const ADT y{2., 2, 1};  // y = 2

  // Test Derivatives with respect to x, i.e. d/dx(Expression)
  EXPECT_FLOAT_EQ((x + y).GetDerivatives()[0], 1.);
  EXPECT_FLOAT_EQ((x * y).GetDerivatives()[0], 2.);
  EXPECT_FLOAT_EQ((x / y).GetDerivatives()[0], 1. / 2.);
  EXPECT_FLOAT_EQ((x - y).GetDerivatives()[0], 1.);

  // Test Derivatives with respect to x, i.e. d/dy(Expression)
  EXPECT_FLOAT_EQ((x + y).GetDerivatives()[1], 1.);
  EXPECT_FLOAT_EQ((x * y).GetDerivatives()[1], 3.);
  EXPECT_FLOAT_EQ((x / y).GetDerivatives()[1], -3. / 4.); // -x/y^2
  EXPECT_FLOAT_EQ((x - y).GetDerivatives()[1], -1.);
};

TEST(AlgoTypeTest, TestValueCorrectnessArithmetic) {
  // Define some Variables
  const ADT x{3., 2, 0};  // x = 3
  const ADT y{2., 2, 1};  // y = 2

  // Test Derivatives with respect to x, i.e. d/dx(Expression)
  EXPECT_FLOAT_EQ(exp(x).GetDerivatives()[0], std::exp(3.));
  EXPECT_FLOAT_EQ(abs(x).GetDerivatives()[0], std::abs(1.));
  EXPECT_FLOAT_EQ(log(x).GetDerivatives()[0], 1. / 3.);
  EXPECT_FLOAT_EQ(log10(x).GetDerivatives()[0], 1. / (3. * std::log(10)));
  EXPECT_FLOAT_EQ(sqrt(x).GetDerivatives()[0], 1. / (2 * std::sqrt(3.)));
  EXPECT_FLOAT_EQ(pow(x, 2.).GetDerivatives()[0], 6.);
  EXPECT_FLOAT_EQ(pow(x, y).GetDerivatives()[0], 2. * 3.);  // y * x^(y-1)
};

TEST(AlgoTypeTest, TestDerivCorrectnessArithmetic) {
  // Define some Variables
  const ADT x{3., 2, 0};  // x = 3
  const ADT y{2., 2, 1};  // y = 2

  EXPECT_FLOAT_EQ(exp(x).GetValue(), std::exp(3.));
  EXPECT_FLOAT_EQ(abs(x).GetValue(), std::abs(3.));
  EXPECT_FLOAT_EQ(log(x).GetValue(), std::log(3.));
  EXPECT_FLOAT_EQ(log10(x).GetValue(), std::log10(3.));
  EXPECT_FLOAT_EQ(sqrt(x).GetValue(), std::sqrt(3.));
  EXPECT_FLOAT_EQ(pow(x, 2.).GetValue(), std::pow(3., 2.));
  EXPECT_FLOAT_EQ(pow(x, y).GetValue(), std::pow(3., 2.));
};

TEST(AlgoTypeTest, TestValueCorrectnessTrigonometric) {
  // Define some Variables
  const ADT x{.7, 2, 0};  // x = .7

  EXPECT_FLOAT_EQ(cos(x).GetValue(), std::cos(.7));
  EXPECT_FLOAT_EQ(sin(x).GetValue(), std::sin(.7));
  EXPECT_FLOAT_EQ(tan(x).GetValue(), std::tan(.7));
  EXPECT_FLOAT_EQ(acos(x).GetValue(), std::acos(.7));
  EXPECT_FLOAT_EQ(asin(x).GetValue(), std::asin(.7));
  EXPECT_FLOAT_EQ(atan(x).GetValue(), std::atan(.7));
};

TEST(AlgoTypeTest, TestDerivCorrectnessTrigonometric) {
  // Define some Variables
  const ADT x{.7, 1, 0};  // x = .7

  EXPECT_FLOAT_EQ(cos(x).GetDerivatives()[0], -std::sin(.7));
  EXPECT_FLOAT_EQ(sin(x).GetDerivatives()[0], std::cos(.7));
  EXPECT_FLOAT_EQ(tan(x).GetDerivatives()[0],
                  1 / (std::cos(.7) * std::cos(.7)));
  EXPECT_FLOAT_EQ(acos(x).GetDerivatives()[0], -1. / std::sqrt(1 - 0.49));
  EXPECT_FLOAT_EQ(asin(x).GetDerivatives()[0], 1. / std::sqrt(1 - 0.49));
  EXPECT_FLOAT_EQ(atan(x).GetDerivatives()[0], 1. / (1. + 0.49));
};

} // namespace beziermanipulation::tests::spline_operations