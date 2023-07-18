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
#include "bezman/src/rational_bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::spline_operations {
using Point3D = Point<3, double>;

class BezierTestingSuite : public ::testing::Test {
 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> line1ctps{Point3D{0., 0., 0.}, Point3D{1., 1., 1.}};
  std::vector<Point3D> lineHOctps{Point3D{0., 0., 0.}, Point3D{1., 1., 0.},
                                  Point3D{2., 0., 0.}};
  std::vector<Point3D> line1actps{Point3D{0., 0., 0.}, Point3D{0.5, 0.5, 0.5}};
  std::vector<Point3D> line1bctps{Point3D{0.5, 0.5, 0.5}, Point3D{1., 1., 1.}};
  std::vector<double> surface_ctps{2., 4., -1., 2.};

  std::array<std::size_t, 1> line_degrees{1};
  std::array<std::size_t, 1> line_ho_degrees{2};
  std::array<std::size_t, 2> surface_degrees{1, 1};

  // Some Lines
  BezierSpline<1, Point3D, double> line1 =
      BezierSpline<1, Point3D, double>(line_degrees, line1ctps);
  BezierSpline<1, Point3D, double> line1a =
      BezierSpline<1, Point3D, double>(line_degrees, line1actps);
  BezierSpline<1, Point3D, double> line1b =
      BezierSpline<1, Point3D, double>(line_degrees, line1bctps);
  BezierSpline<2, double, double> surface{surface_degrees, surface_ctps};
  BezierSpline<1, Point3D, double> lineHO{line_ho_degrees, lineHOctps};
};

// Compare two known splines
TEST_F(BezierTestingSuite, TestSplittingLine) {
  // Expect equality.
  EXPECT_EQ(line1.SplitAtPosition(0.5)[0], line1a);
  EXPECT_EQ(line1.SplitAtPosition(0.5)[1], line1b);
}

// Compare two lines at specific points with higher order
TEST_F(BezierTestingSuite, TestSplittingLineHigherOrder) {
  const auto split_position = 0.5;
  const auto splitline = lineHO.SplitAtPosition(split_position)[0];
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(splitline.Evaluate(x)[0],
                    lineHO.Evaluate(x * split_position)[0]);
  }
}

// Compare distinct samples
TEST_F(BezierTestingSuite, TestSplittingSurface) {
  // Expect equality.
  const double split_plane =
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  const auto surface_split_first_half =
      surface.SplitAtPosition(split_plane, 1)[0];
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(surface_split_first_half.Evaluate(x, y),
                    surface.Evaluate(x, split_plane * y));
  }
}

// Compare distinct samples
TEST_F(BezierTestingSuite, TestSplittingSurfaceHighOrder) {
  // Expect equality.
  const double split_plane =
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  surface.OrderElevateAlongParametricDimension(0)
      .OrderElevateAlongParametricDimension(1);
  const auto surface_split_first_half =
      surface.SplitAtPosition(split_plane, 1)[0];
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(surface_split_first_half.Evaluate(x, y),
                    surface.Evaluate(x, split_plane * y));
  }
}

}  // namespace bezman::tests::spline_operations