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

#include "bezman/src/utils/algorithms/bernstein_polynomial.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace bezman;

namespace bezman::tests::bernstein_polynomial_test {

class BernsteinPolynomialTest : public ::testing::Test {
 protected:
  void SetUp() override {}

  const double PI = std::acos(-1);

  // Analytical implementation of some polynomials
  double b02(const double& x) { return 1 - 2 * x + x * x; }
  double b12(const double& x) { return 2 * x - 2 * x * x; }
  double b22(const double& x) { return x * x; }
  double b02_1(const double& x) { return -2. + 2. * x; }
  double b12_1(const double& x) { return 2. - 4. * x; }
  double b22_1(const double& x) { return 2. * x; }
  double b02_2([[maybe_unused]] const double& x) { return 2.; }
  double b12_2([[maybe_unused]] const double& x) { return -4.; }
  double b22_2([[maybe_unused]] const double& x) { return 2.; }

  // Envelop function (see )
  double BernsteinEnvelope(const std::size_t& degree, const double& value) {
    return 1. / std::sqrt(2 * PI * degree * value * (1 - value));
  }
};

TEST_F(BernsteinPolynomialTest, TestEnvelopAndPartitionOfUnity) {
  std::size_t n_tests = 10;
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    const std::size_t degree = rand() % 30;
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto evaluations =
        bezman::utils::algorithms::BernsteinPolynomial<double>::Evaluate(degree,
                                                                         x);
    double sum{};
    const double env_v = BernsteinEnvelope(degree, x);
    for (std::size_t i_basis{}; i_basis < degree + 1; i_basis++) {
      sum += evaluations[i_basis];
      EXPECT_TRUE(evaluations[i_basis] <= env_v);
    }
    EXPECT_FLOAT_EQ(sum, 1.);
  }
}

TEST_F(BernsteinPolynomialTest, TestDerivatives) {
  std::size_t n_tests = 10;
  std::size_t max_degree = 5;
  // Test that derivatives sum up to 0
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    const std::size_t degree = rand() % max_degree + 2;
    // Make sure deriv is at least 1 but <= degree
    const std::size_t deriv = rand() % (degree - 1) + 1;
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto deriv_evaluations =
        bezman::utils::algorithms::BernsteinPolynomial<
            double>::EvaluateDerivative(degree, x, deriv);
    double sum{};
    for (std::size_t i_basis{}; i_basis < (degree + 1); i_basis++) {
      sum += deriv_evaluations[i_basis];
    }
    // The sum should be zero, however, the values are extremely high for high
    // derivatives, so it should be scaled (order of magnitude is n!/(n-p!),)
    EXPECT_TRUE(sum < 1e-10);
  }

  // Test explicitely for second degree polynomial
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    // Make sure deriv is at least 1 but <= degree
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto evaluations = bezman::utils::algorithms::BernsteinPolynomial<
        double>::EvaluateDerivative(2, x, 0);
    const auto deriv_evaluations =
        bezman::utils::algorithms::BernsteinPolynomial<
            double>::EvaluateDerivative(2, x, 1);
    const auto sec_deriv_evaluations =
        bezman::utils::algorithms::BernsteinPolynomial<
            double>::EvaluateDerivative(2, x, 2);
    EXPECT_FLOAT_EQ(evaluations[0], b02(x));
    EXPECT_FLOAT_EQ(evaluations[1], b12(x));
    EXPECT_FLOAT_EQ(evaluations[2], b22(x));
    EXPECT_FLOAT_EQ(deriv_evaluations[0], b02_1(x));
    EXPECT_FLOAT_EQ(deriv_evaluations[1], b12_1(x));
    EXPECT_FLOAT_EQ(deriv_evaluations[2], b22_1(x));
    EXPECT_FLOAT_EQ(sec_deriv_evaluations[0], b02_2(x));
    EXPECT_FLOAT_EQ(sec_deriv_evaluations[1], b12_2(x));
    EXPECT_FLOAT_EQ(sec_deriv_evaluations[2], b22_2(x));
  }
}

}  // namespace bezman::tests::bernstein_polynomial_test
