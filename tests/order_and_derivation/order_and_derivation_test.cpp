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
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point2D> line1ctps{Point2D{0., 0.}, Point2D{1., 1.}};
  std::vector<Point2D> line1_order_elev_ctps{Point2D{0., 0.}, Point2D{0.5, 0.5},
                                             Point2D{1., 1.}};
  std::vector<Point2D> line1_derv_ctps{Point2D{1., 1.}};
  std::vector<Point2D> surface_ctps{
      Point2D{0., 0.},    Point2D{0.5, 0.2}, Point2D{1., 0.},
      Point2D{-0.2, 0.5}, Point2D{0.5, 0.5}, Point2D{1.2, 0.5},
      Point2D{0., 1.},    Point2D{0.5, 0.8}, Point2D{1., 1.},
  };

  std::array<std::size_t, 1> degrees_line{1};
  std::array<std::size_t, 1> degrees_line_order_elev{2};
  std::array<std::size_t, 1> degrees_deriv{0};
  std::array<std::size_t, 2> surface_degrees{2, 2};

  // Some Lines
  BezierSpline<1, Point2D, double> line1 =
      BezierSpline<1, Point2D, double>(degrees_line, line1ctps);
  BezierSpline<1, Point2D, double> line1_deriv =
      BezierSpline<1, Point2D, double>(degrees_deriv, line1_derv_ctps);
  BezierSpline<1, Point2D, double> line1_order_eliv =
      BezierSpline<1, Point2D, double>(degrees_line_order_elev,
                                       line1_order_elev_ctps);

  // Some second order scewed plane
  BezierSpline<2, Point2D, double> surface_spline{surface_degrees,
                                                  surface_ctps};

  auto CreateRandomSpline(unsigned int degree) {
    BezierSpline<1, double, double> randomSpline{
        std::array<std::size_t, 1>{degree}};
    for (unsigned int i{}; i < degree; i++) {
      randomSpline.ControlPoint(i) =
          static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }
    return randomSpline;
  }
};

// Order Elevation Test
TEST_F(BezierTestingSuite, OrderElevation) {
  // Expect equality.
  EXPECT_EQ(line1.OrderElevateAlongParametricDimension(0), line1_order_eliv);
}

// Demonstrate Elevation at random Points and random lines
TEST_F(BezierTestingSuite, OrderElevation2) {
  // Expect equality.
  auto spline1 = CreateRandomSpline(10);
  const auto spline1_original = spline1;
  // Elevate the order a couple of times
  spline1.OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(0);

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.Evaluate(x), spline1_original.Evaluate(x));
  }
}

// Demonstrate Elevation at random Points and random lines
TEST_F(BezierTestingSuite, OrderElevation3) {
  // Expect equality.
  const auto surface_copy{surface_spline};
  // Elevate the order a couple of times
  surface_spline.OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(1)
      .OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(1)
      .OrderElevateAlongParametricDimension(1);

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const double y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(surface_copy.Evaluate(x, y)[0],
                    surface_spline.Evaluate(x, y)[0]);
    EXPECT_FLOAT_EQ(surface_copy.Evaluate(x, y)[1],
                    surface_spline.Evaluate(x, y)[1]);
  }
}

// Derivation Of line with known results
TEST_F(BezierTestingSuite, Derivation) {
  // Expect equality.
  EXPECT_EQ(line1.DerivativeWRTParametricDimension(0), line1_deriv);
}
}  // namespace bezman::tests::basic_operations