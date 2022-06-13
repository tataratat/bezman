#include <gtest/gtest.h>

#include <array>

#include "bezman/src/bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::basic_operations {

class BezierTestingSuite : public ::testing::Test {
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

  const auto CreateRandomSpline(unsigned int degree) {
    BezierSpline<1, double, double> randomSpline{
        std::array<std::size_t, 1>{degree}};
    for (unsigned int i{}; i < degree; i++) {
      randomSpline.ControlPoint(i) =
          static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
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

}  // namespace bezman::tests::basic_operations