#include <gtest/gtest.h>

#include <array>

#define MAX_BINOMIAL_DEGREE 62u // Required for second test set
#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::composition_test {

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
  const auto composed_spline = surface.compose(line);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto composed_spline_position = composed_spline.evaluate(x);
    const auto composed_spline_expectation = surface.evaluate(line.evaluate(x));
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
    surface.order_elevate_along_parametric_dimension(0)
        .order_elevate_along_parametric_dimension(1);
    line.order_elevate_along_parametric_dimension(0);
  }
  const auto composed_spline = surface.compose(line);

  std::cerr << "[          ] surface-degrees: =     \t(" << surface.degrees[0]
            << ", " << surface.degrees[1] << ")" << std::endl;
  std::cerr << "[          ] line-degree: =         \t(" << line.degrees[0]
            << ")" << std::endl;

  std::cerr << "[          ] composition-degrees: = \t("
            << composed_spline.degrees[0] << ")" << std::endl;
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto composed_spline_position = composed_spline.evaluate(x);
    const auto composed_spline_expectation = surface.evaluate(line.evaluate(x));
    for (unsigned int i_dim{}; i_dim < 3; i_dim++) {
      // Compare the different dimensions element-wise
      EXPECT_FLOAT_EQ(composed_spline_position[i_dim],
                      composed_spline_expectation[i_dim]);
    }
  }
}

}  // namespace beziermanipulation::tests::composition_test