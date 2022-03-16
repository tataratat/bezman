#include <gtest/gtest.h>

#include <array>

#include "bezierManipulation/src/bezier_spline_group.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::spline_operations {

class BezierTestingSuite : public ::testing::Test {
  using Point3D = Point<3, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> line1ctps{Point3D{0., 0., 0.}, Point3D{1., 1., 1.}};
  std::vector<Point3D> line1actps{Point3D{0., 0., 0.}, Point3D{0.5, 0.5, 0.5}};
  std::vector<Point3D> line1bctps{Point3D{0.5, 0.5, 0.5}, Point3D{1., 1., 1.}};
  std::vector<double> surface_ctps{2., 4., -1., 2.};

  std::array<std::size_t, 1> line_degrees{1};
  std::array<std::size_t, 2> surface_degrees{1, 1};

  // Some Lines
  BezierSpline<1, Point3D, double> line1 =
      BezierSpline<1, Point3D, double>(line_degrees, line1ctps);
  BezierSpline<1, Point3D, double> line1a =
      BezierSpline<1, Point3D, double>(line_degrees, line1actps);
  BezierSpline<1, Point3D, double> line1b =
      BezierSpline<1, Point3D, double>(line_degrees, line1bctps);
  BezierSpline<2, double, double> surface =
      BezierSpline<2, double, double>(surface_degrees, surface_ctps);
};

// Compare two known splines
TEST_F(BezierTestingSuite, TestSplittingLine) {
  // Expect equality.
  EXPECT_EQ(line1.split(0.5)[0], line1a);
  EXPECT_EQ(line1.split(0.5)[1], line1b);
}

// Compare distinct samples
TEST_F(BezierTestingSuite, TestSplittingSurface) {
  // Expect equality.
  const double split_plane =
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  const auto surface_split_first_half = surface.split(split_plane, 1)[0];
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(surface_split_first_half.evaluate(x, y),
                    surface.evaluate(x, split_plane * y));
  }
};
}  // namespace beziermanipulation::tests::spline_operations