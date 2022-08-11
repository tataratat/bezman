/*
MIT License

Copyright (c) 2022 zwar@ilsb.tuwien.ac.at

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "bezman/src/utils/computational_differentiation/algo_diff_type.hpp"

#include <gtest/gtest.h>

#include <cmath>

using namespace bezman;

namespace bezman::tests::algo_diff_type_test {

using ADT = utils::computational_differentiation::AlgoDiffType<double>;
using ADTstatic = utils::computational_differentiation::AlgoDiffType<double, 2>;

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
}

// Constructors and Basic Operations
TEST(AlgoTypeTest, TestValueCorrectnessBasicFriends) {
  // Define some Variables
  const double x{3.};  // x = 3
  ADT y{2., 1, 0};     // y = 2

  // Test Values
  EXPECT_FLOAT_EQ((x + y).GetValue(), 3. + 2.);
  EXPECT_FLOAT_EQ((x * y).GetValue(), 3. * 2.);
  EXPECT_FLOAT_EQ((x / y).GetValue(), 3. / 2.);
  EXPECT_FLOAT_EQ((x - y).GetValue(), 3. - 2.);
}

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
  EXPECT_FLOAT_EQ((x / y).GetDerivatives()[1], -3. / 4.);  // -x/y^2
  EXPECT_FLOAT_EQ((x - y).GetDerivatives()[1], -1.);
}

TEST(AlgoTypeTest, TestDerivCorrectnessBasicFriends) {
  // Define some Variables
  double x{3.};           // x = 3
  const ADT y{2., 1, 0};  // y = 2

  // Test Derivatives with respect to x, i.e. d/dx(Expression)
  EXPECT_FLOAT_EQ((x + y).GetDerivatives()[0], 1.);
  EXPECT_FLOAT_EQ((x * y).GetDerivatives()[0], 3.);
  EXPECT_FLOAT_EQ((x / y).GetDerivatives()[0], -3. / 4.);
  EXPECT_FLOAT_EQ((x - y).GetDerivatives()[0], -1.);
}

TEST(AlgoTypeTest, TestDerivCorrectnessArithmetic) {
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
}

TEST(AlgoTypeTest, TestValueCorrectnessArithmetic) {
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
}

TEST(AlgoTypeTest, TestValueCorrectnessTrigonometric) {
  // Define some Variables
  const ADT x{.7, 2, 0};  // x = .7

  EXPECT_FLOAT_EQ(cos(x).GetValue(), std::cos(.7));
  EXPECT_FLOAT_EQ(sin(x).GetValue(), std::sin(.7));
  EXPECT_FLOAT_EQ(tan(x).GetValue(), std::tan(.7));
  EXPECT_FLOAT_EQ(acos(x).GetValue(), std::acos(.7));
  EXPECT_FLOAT_EQ(asin(x).GetValue(), std::asin(.7));
  EXPECT_FLOAT_EQ(atan(x).GetValue(), std::atan(.7));
}

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
}

// Constructors and Basic Operations for static types
TEST(AlgoTypeTest, TestValueCorrectnessBasicStatic) {
  // Define some Variables
  ADTstatic x{3., 0};  // x = 3
  ADTstatic y{2., 1};  // y = 2

  // Test Values
  EXPECT_FLOAT_EQ((x + y).GetValue(), 3. + 2.);
  EXPECT_FLOAT_EQ((x * y).GetValue(), 3. * 2.);
  EXPECT_FLOAT_EQ((x / y).GetValue(), 3. / 2.);
  EXPECT_FLOAT_EQ((x - y).GetValue(), 3. - 2.);
}

TEST(AlgoTypeTest, TestDerivCorrectnessTrigonometricStatic) {
  // Define some Variables
  const ADTstatic x{.7, 0};  // x = .7

  EXPECT_FLOAT_EQ(cos(x).GetDerivatives()[0], -std::sin(.7));
  EXPECT_FLOAT_EQ(sin(x).GetDerivatives()[0], std::cos(.7));
  EXPECT_FLOAT_EQ(tan(x).GetDerivatives()[0],
                  1 / (std::cos(.7) * std::cos(.7)));
  EXPECT_FLOAT_EQ(acos(x).GetDerivatives()[0], -1. / std::sqrt(1 - 0.49));
  EXPECT_FLOAT_EQ(asin(x).GetDerivatives()[0], 1. / std::sqrt(1 - 0.49));
  EXPECT_FLOAT_EQ(atan(x).GetDerivatives()[0], 1. / (1. + 0.49));
}

}  // namespace bezman::tests::algo_diff_type_test