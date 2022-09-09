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

#include "bezman/src/bezier_group.hpp"

#include <gtest/gtest.h>

#include <array>

#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/rational_bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::bezier_group_test {

class BezierGroupTest : public ::testing::Test {
 public:
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;
  using BezierLine3D = BezierSpline<1, Point3D, double>;
  using RationalBezierLine3D = RationalBezierSpline<1, Point3D, double>;

 protected:
  void SetUp() override {}

  std::vector<Point3D> line1ctps{Point3D{0., 0., 0.}, Point3D{0., 0.5, 0.},
                                 Point3D{0., 1., 0.}};
  std::array<std::size_t, 1> degrees{2};

  BezierLine3D bezier_line_spline{degrees, line1ctps};
  BezierGroup<BezierLine3D> bezier_line_group{bezier_line_spline};

  RationalBezierLine3D rational_line_spline{bezier_line_spline};
  BezierGroup<RationalBezierLine3D> rational_line_group{rational_line_spline};

  // Higher Dimensional Examples
  std::vector<double> weights_arc{1., 1 / std::sqrt(2), 1.};
  std::vector<Point2D> ctps_arc{Point2D{1., 0.}, Point2D{1., 1.},
                                Point2D{0., 1.}};

  // Create Surfaces
  std::array<std::size_t, 2> degrees_2D{2, 1};
  std::vector<Point2D> ctps_2D{
      Point2D{1., 0.}, Point2D{1., 1.}, Point2D{0., 1.},
      Point2D{2., 0.}, Point2D{2., 2.}, Point2D{0., 2.},
  };
  std::vector<double> weights_2D{1., 1 / std::sqrt(2), 1.,
                                 1., 1 / std::sqrt(2), 1.};
};

TEST_F(BezierGroupTest, TestFitIntoUnitCubeFunction) {
  EXPECT_TRUE(bezier_line_group.FitsIntoUnitCube());
  EXPECT_TRUE(rational_line_group.FitsIntoUnitCube());

  // Create new spline that is too big and fit it
  auto bezier_too_big = bezier_line_spline * 2.;
  EXPECT_NO_FATAL_FAILURE(bezier_line_group.push_back(bezier_too_big));
  EXPECT_FALSE(bezier_line_group.FitsIntoUnitCube());
  EXPECT_NO_FATAL_FAILURE(bezier_line_group.FitIntoUnitCube());
  EXPECT_TRUE(bezier_line_group.FitsIntoUnitCube());
}

TEST_F(BezierGroupTest, TestMaximumMinimumCornerFunction) {
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  const RationalBezierSurface circle_segment{degrees_2D, ctps_2D, weights_2D};
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degrees, ctps_arc, weights_arc);
  const BezierGroup<RationalBezier> circular_arc_group =
      BezierGroup<RationalBezier>(circular_arc);

  EXPECT_FLOAT_EQ(circle_segment.MaximumCorner()[0], Point2D(2., 2.)[0]);
  EXPECT_FLOAT_EQ(circle_segment.MaximumCorner()[1], Point2D(2., 2.)[1]);
  EXPECT_FLOAT_EQ(circle_segment.MinimumCorner()[0], Point2D(0., 0.)[0]);
  EXPECT_FLOAT_EQ(circle_segment.MinimumCorner()[1], Point2D(0., 0.)[1]);
}

TEST_F(BezierGroupTest, TestComposeSplineGroupFunction) {
  using RationalBezierLine = RationalBezierSpline<1, Point2D, double>;
  using BezierLine = BezierSpline<1, Point2D, double>;
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  using BezierSurface = BezierSpline<2, Point2D, double>;

  const BezierLine bezierLine{degrees, ctps_arc};
  const RationalBezierLine rational_line(bezierLine);
  const RationalBezierSurface circle_segment_rat{degrees_2D, ctps_2D,
                                                 weights_2D};
  const BezierSurface circle_segment_pol{degrees_2D, ctps_2D};

  BezierGroup<RationalBezierSurface> surface_rat{circle_segment_rat};
  BezierGroup<BezierSurface> surface_pol{circle_segment_pol};
  BezierGroup<RationalBezierLine> line_rat{rational_line};
  BezierGroup<BezierLine> line_pol{bezierLine};

  EXPECT_NO_FATAL_FAILURE(surface_pol.Compose(line_pol));
  EXPECT_NO_FATAL_FAILURE(surface_rat.Compose(line_rat));
  EXPECT_NO_FATAL_FAILURE(surface_pol.Compose(line_rat));
  EXPECT_NO_FATAL_FAILURE(surface_rat.Compose(line_pol));

  // Test Results
  for (int i{}; i < 2; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    for (unsigned int i_dim{}; i_dim < 2; i_dim++) {
      EXPECT_FLOAT_EQ(surface_pol.Compose(line_pol)[0].Evaluate(x)[i_dim],
                      surface_pol[0].Evaluate(line_pol[0].Evaluate(x))[i_dim]);
      EXPECT_FLOAT_EQ(surface_rat.Compose(line_rat)[0].Evaluate(x)[i_dim],
                      surface_rat[0].Evaluate(line_rat[0].Evaluate(x))[i_dim]);
      EXPECT_FLOAT_EQ(surface_pol.Compose(line_rat)[0].Evaluate(x)[i_dim],
                      surface_pol[0].Evaluate(line_rat[0].Evaluate(x))[i_dim]);
      EXPECT_FLOAT_EQ(surface_rat.Compose(line_pol)[0].Evaluate(x)[i_dim],
                      surface_rat[0].Evaluate(line_pol[0].Evaluate(x))[i_dim]);
    }
  }
}

}  // namespace bezman::tests::bezier_group_test