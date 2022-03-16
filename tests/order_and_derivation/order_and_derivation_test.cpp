#include <gtest/gtest.h>

#include <array>

#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::basic_operations {

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
  std::array<std::size_t, 1> degrees_line{1};
  std::array<std::size_t, 1> degrees_line_order_elev{2};
  std::array<std::size_t, 1> degrees_deriv{0};

  // Some Lines
  BezierSpline<1, Point2D, double> line1 =
      BezierSpline<1, Point2D, double>(degrees_line, line1ctps);
  BezierSpline<1, Point2D, double> line1_deriv =
      BezierSpline<1, Point2D, double>(degrees_deriv, line1_derv_ctps);
  BezierSpline<1, Point2D, double> line1_order_eliv =
      BezierSpline<1, Point2D, double>(degrees_line_order_elev,
                                       line1_order_elev_ctps);

  auto CreateRandomSpline(unsigned int degree){
    BezierSpline<1, double, double> randomSpline{std::array<std::size_t, 1>{degree}};
    for (unsigned int i{}; i < degree; i++){
      randomSpline.control_point(i) =
          static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }
    return randomSpline;
  }
};

// Order Elevation Test
TEST_F(BezierTestingSuite, OrderElevation) {
  // Expect equality.
  EXPECT_EQ(line1.order_elevate_along_parametric_dimension(0), line1_order_eliv);
}

// Demonstrate Elevation at random Points and random lines
TEST_F(BezierTestingSuite, OrderElevation2) {
  // Expect equality.
  auto spline1 = CreateRandomSpline(10);
  const auto spline1_original = spline1;
  // Elevate the order a couple of times
  spline1.order_elevate_along_parametric_dimension(0)
      .order_elevate_along_parametric_dimension(0)
      .order_elevate_along_parametric_dimension(0)
      .order_elevate_along_parametric_dimension(0)
      .order_elevate_along_parametric_dimension(0);

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.evaluate(x), spline1_original.evaluate(x));
  }
}

// Derivation Of line with known results
TEST_F(BezierTestingSuite, Derivation) {
  // Expect equality.
  EXPECT_EQ(line1.derive_along_parametric_dimension(0), line1_deriv);
}
}  // namespace beziermanipulation::tests::basic_operations