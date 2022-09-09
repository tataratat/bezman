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

#define MAX_BINOMIAL_DEGREE 62u  // Required for second test set
#include "bezman/src/bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::composition_test {

class BezierTestingSuite : public ::testing::Test {
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> surface_ctps{Point3D{0., 0., 0.}, Point3D{2., 1., 0.},
                                    Point3D{0., 2., 0.}, Point3D{2., 3., 0.}};
  std::vector<Point2D> line_ctps{Point2D{0., 0.}, Point2D{1., 0.},
                                 Point2D{1., 1.}};
  std::array<std::size_t, 2> surface_degrees{1, 1};
  std::array<std::size_t, 1> line_degrees{2};

  // Define Splines
  BezierSpline<2, Point3D, double> surface =
      BezierSpline<2, Point3D, double>(surface_degrees, surface_ctps);
  BezierSpline<1, Point2D, double> line =
      BezierSpline<1, Point2D, double>(line_degrees, line_ctps);
};

/*
 * Demonstrate spline composition
 *
 * There is no analytical example implemented, where the spline ctps are
 * explicitly known. Thus, the composed spline is tested at random points within
 * the domain
 */
TEST_F(BezierTestingSuite, Composition) {
  // Expect equality.
  const auto composed_spline = surface.Compose(line);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto composed_spline_position = composed_spline.Evaluate(x);
    const auto composed_spline_expectation = surface.Evaluate(line.Evaluate(x));
    for (unsigned int i_dim{}; i_dim < 3; i_dim++) {
      // Compare the different dimensions element-wise
      EXPECT_FLOAT_EQ(composed_spline_position[i_dim],
                      composed_spline_expectation[i_dim]);
    }
  }
}

/*
 * Demonstrate some basic compositions on a more expensive example with the same
 * geometry
 */
TEST_F(BezierTestingSuite, CompositionHighOrder) {
  // Increase the order a few times
  for (int i{}; i < 4; i++) {
    surface.OrderElevateAlongParametricDimension(0)
        .OrderElevateAlongParametricDimension(1);
    line.OrderElevateAlongParametricDimension(0);
  }
  const auto composed_spline = surface.Compose(line);

  std::cerr << "[          ] surface-degrees: =     \t("
            << surface.GetDegrees()[0] << ", " << surface.GetDegrees()[1] << ")"
            << std::endl;
  std::cerr << "[          ] line-degree: =         \t(" << line.GetDegrees()[0]
            << ")" << std::endl;

  std::cerr << "[          ] composition-degrees: = \t("
            << composed_spline.GetDegrees()[0] << ")" << std::endl;
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto composed_spline_position = composed_spline.Evaluate(x);
    const auto composed_spline_expectation = surface.Evaluate(line.Evaluate(x));
    for (unsigned int i_dim{}; i_dim < 3; i_dim++) {
      // Compare the different dimensions element-wise
      EXPECT_FLOAT_EQ(composed_spline_position[i_dim],
                      composed_spline_expectation[i_dim]);
    }
  }
}

}  // namespace bezman::tests::composition_test