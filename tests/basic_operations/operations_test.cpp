#include <gtest/gtest.h>

#include <array>

#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::basic_operations {

class BezierBasicOperationSuite : public ::testing::Test {
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

  const auto CreateRandomSpline(unsigned int degree){
    BezierSpline<1, double, double> randomSpline{std::array<std::size_t, 1>{degree}};
    for (unsigned int i{}; i < degree; i++){
      randomSpline.control_point(i) =
          static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }
    return randomSpline;
  }
};

// Demonstrate some basic assertions.
TEST_F(BezierBasicOperationSuite, TestAddition) {
  // Expect equality.
  EXPECT_EQ(line1 + line2, line3);
  EXPECT_EQ(line1_copy.order_elevate_along_parametric_dimension(0) + line2,
            line3_copy.order_elevate_along_parametric_dimension(0));
  EXPECT_EQ(line1 + line2_copy.order_elevate_along_parametric_dimension(0),
            line3.order_elevate_along_parametric_dimension(0));
}

// Demonstrate some basic assertions.
TEST_F(BezierBasicOperationSuite, TestAddition2) {
  // Expect equality.

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_EQ(line1.evaluate(x) + line2.evaluate(x),
              (line1 + line2).evaluate(x));
  }
}

// Demonstrate some basic assertions.
TEST_F(BezierBasicOperationSuite, TestAddition3) {
  // Expect equality.
  const auto spline1 = CreateRandomSpline(10);
  const auto spline2 = CreateRandomSpline(15);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.evaluate(x) + spline2.evaluate(x),
                    (spline1 + spline2).evaluate(x));
  }
}

// Demonstrate some basic assertions.
TEST_F(BezierBasicOperationSuite, MultiplicationTest1) {
  // Expect equality.

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(line1.evaluate(x) * line2.evaluate(x),
                    (line1 * line2).evaluate(x));
  }
}

// Demonstrate some basic assertions.
TEST_F(BezierBasicOperationSuite, MultiplicationTest2) {
  // Expect equality.
  const auto spline1 = CreateRandomSpline(10);
  const auto spline2 = CreateRandomSpline(15);
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(spline1.evaluate(x) * spline2.evaluate(x),
                    (spline1 * spline2).evaluate(x));
  }
}

}  // namespace beziermanipulation::tests::basic_operations