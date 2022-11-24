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

#include <gtest/gtest.h>

#include <array>

#include "bezman/src/bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::basic_operations {

class BezierTestingSuite : public ::testing::Test {
 public:
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> line1ctps{Point3D{0., 0., 0.}, Point3D{0., 0.5, 0.},
                                 Point3D{0., 1., 0.}};
  std::vector<Point3D> line2ctps{Point3D{0., 0., 0.}, Point3D{0.5, 0., 0.},
                                 Point3D{1., 0., 0.}};
  std::vector<Point3D> line3ctps{Point3D{0., 0., 0.}, Point3D{0.5, 0.5, 0.},
                                 Point3D{1., 1., 0.}};
  std::array<std::size_t, 1> degrees{2};

  // Some Lines
  BezierSpline<1, Point3D, double> line1 =
      BezierSpline<1, Point3D, double>(degrees, line1ctps);
  BezierSpline<1, Point3D, double> line2 =
      BezierSpline<1, Point3D, double>(degrees, line2ctps);
  BezierSpline<1, Point3D, double> line3 =
      BezierSpline<1, Point3D, double>(degrees, line3ctps);
  BezierSpline<1, Point3D, double> line1_copy{line1};
  BezierSpline<1, Point3D, double> line2_copy{line2};
  BezierSpline<1, Point3D, double> line3_copy{line3};

  // Surface for testing evaluation routines
  std::vector<Point2D> surface_ctps{
      Point2D{0., 0.},    Point2D{0.5, 0.2}, Point2D{1., 0.},
      Point2D{-0.2, 0.5}, Point2D{0.5, 0.5}, Point2D{1.2, 0.5},
      Point2D{0., 1.},    Point2D{0.5, 0.8}, Point2D{1., 1.},
  };
  std::array<std::size_t, 2> surface_degrees{2, 2};
  BezierSpline<2, Point2D, double> surface_spline{surface_degrees,
                                                  surface_ctps};
  template <typename PhysicalPoint = double>
  auto CreateRandomSpline(unsigned int degree) {
    BezierSpline<1, PhysicalPoint, double> randomSpline{
        std::array<std::size_t, 1>{degree}};

    for (unsigned int i{}; i < degree; i++) {
      if constexpr (std::is_scalar_v<PhysicalPoint>) {
        randomSpline.ControlPoint(i) =
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
      } else {
        for (std::size_t i_dim{}; i_dim < PhysicalPoint::kSpatialDimension;
             i_dim++) {
          randomSpline.ControlPoint(i)[i_dim] =
              static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        }
      }
    }
    return randomSpline;
  }
};

// Demonstrate Additions with known results
TEST_F(BezierTestingSuite, TestAddition) {
  // Expect equality.
  EXPECT_EQ(line1 + line2, line3);
  EXPECT_EQ(line1_copy.OrderElevateAlongParametricDimension(0) + line2,
            line3_copy.OrderElevateAlongParametricDimension(0));
  EXPECT_EQ(line1 + line2_copy.OrderElevateAlongParametricDimension(0),
            line3.OrderElevateAlongParametricDimension(0));
  EXPECT_EQ(-line1, (-1) * line1);
}

// Demonstrate Additions at random points
TEST_F(BezierTestingSuite, TestAddition2) {
  // Expect equality.

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ((line1.Evaluate(x) + line2.Evaluate(x))[0],
                    (line1 + line2).ForwardEvaluate(x)[0]);
    EXPECT_FLOAT_EQ((line1.Evaluate(x) + line2.Evaluate(x))[1],
                    (line1 + line2).ForwardEvaluate(x)[1]);
    EXPECT_FLOAT_EQ((line1.Evaluate(x) + line2.Evaluate(x))[2],
                    (line1 + line2).ForwardEvaluate(x)[2]);
  }
}

// Demonstrate Additions at random points
TEST_F(BezierTestingSuite, TestEvaluationRoutines) {
  // Expect equalityS.
  for (int i{0}; i < 9; i++) {
    surface_spline.OrderElevateAlongParametricDimension(0);
    surface_spline.OrderElevateAlongParametricDimension(1);
  }

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const double y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(surface_spline.Evaluate(x, y)[0],
                    surface_spline.ForwardEvaluate(x, y)[0]);
    EXPECT_FLOAT_EQ(surface_spline.Evaluate(x, y)[1],
                    surface_spline.ForwardEvaluate(x, y)[1]);
  }
}

// Test Basis Function evaluation
TEST_F(BezierTestingSuite, TestBasisFunctionContributions) {
  // Elevate dergees
  surface_spline.OrderElevateAlongParametricDimension(0);
  surface_spline.OrderElevateAlongParametricDimension(1);
  const auto degrees = surface_spline.GetDegrees();
  // Create random evaluation point
  const double xx{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  const double xy{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  const Point2D x{xx, xy};
  // Retrieve basis functions
  const auto basis_functions =
      surface_spline.BasisFunctionContributions(xx, xy);

  // Compare to analytical solution
  for (std::size_t pdim{}; pdim < 2; pdim++) {
    for (std::size_t i_basis{}; i_basis < degrees[pdim] + 1; i_basis++) {
      const double analytical_solution_bernstein_pol =
          utils::FastBinomialCoefficient::choose(degrees[pdim], i_basis) *
          std::pow(x[pdim], i_basis) *
          std::pow(1. - x[pdim], degrees[pdim] - i_basis);
      EXPECT_FLOAT_EQ(basis_functions[pdim][i_basis],
                      analytical_solution_bernstein_pol);
    }
  }
}

// Test Basis Function Evaluations
TEST_F(BezierTestingSuite, TestBasisFunctions) {
  // Elevate dergees
  surface_spline.OrderElevateAlongParametricDimension(0);
  surface_spline.OrderElevateAlongParametricDimension(1);
  const auto degrees = surface_spline.GetDegrees();
  const std::size_t n_tests{5};
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    // Create random evaluation point
    const double xx{static_cast<double>(rand()) /
                    static_cast<double>(RAND_MAX)};
    const double xy{static_cast<double>(rand()) /
                    static_cast<double>(RAND_MAX)};
    // Retrieve basis functions
    const auto basis_functions = surface_spline.BasisFunctions(xx, xy);

    // Compare to analytical solution
    for (std::size_t i_basis{}; i_basis < degrees[0] + 1; i_basis++) {
      for (std::size_t j_basis{}; j_basis < degrees[1] + 1; j_basis++) {
        const double analytical_solution_bernstein_pol =
            utils::FastBinomialCoefficient::choose(degrees[0], i_basis) *
            std::pow(xx, i_basis) * std::pow(1. - xx, degrees[0] - i_basis) *
            utils::FastBinomialCoefficient::choose(degrees[1], j_basis) *
            std::pow(xy, j_basis) * std::pow(1. - xy, degrees[1] - j_basis);
        EXPECT_FLOAT_EQ(basis_functions[i_basis + (degrees[0] + 1) * j_basis],
                        analytical_solution_bernstein_pol);
      }
    }
  }
}

// Demonstrate Addition at random Points and random lines
TEST_F(BezierTestingSuite, TestAddition3) {
  // Expect equality.
  const auto spline1 = CreateRandomSpline(10);
  const auto spline2 = CreateRandomSpline(15);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.Evaluate(x) + spline2.Evaluate(x),
                    (spline1 + spline2).Evaluate(x));
  }
}

// Multiplications at random points
TEST_F(BezierTestingSuite, MultiplicationTest1) {
  // Expect equality.

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(line1.Evaluate(x) * line2.Evaluate(x),
                    (line1 * line2).Evaluate(x));
  }
}
// Multiplications at random points and random splines
TEST_F(BezierTestingSuite, MultiplicationTest2) {
  // Expect equality.
  const auto spline1 = CreateRandomSpline(10);
  const auto spline2 = CreateRandomSpline(15);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.Evaluate(x) * spline2.Evaluate(x),
                    (spline1 * spline2).Evaluate(x));
  }
}

// Raise a random spline to some powers and test the results at random points to
// check for equality
TEST_F(BezierTestingSuite, RaisingPowersTest) {
  // Expect equality.
  const auto spline_base = CreateRandomSpline(2);
  const std::size_t maximum_power_test{5};
  for (std::size_t i{}; i < maximum_power_test; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ((spline_base.RaisePower(i)).Evaluate(x),
                    std::pow(spline_base.Evaluate(x), i));
  }
}

}  // namespace bezman::tests::basic_operations