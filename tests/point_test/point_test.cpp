#include <gtest/gtest.h>

#include <array>

#define MAX_BINOMIAL_DEGREE 62u  // Required for second test set
#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/utils/type_traits/is_point.hpp"

using namespace bezman;

namespace bezman::tests::point_test {

class BezierTestingSuite : public ::testing::Test {
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Using ints to ensure precise comparison.
  Point<3, int> point1{1, 2, 3};
  Point<3, int> point2{2, 3, 2};
  Point<3, int> point3{-1, -1, 1};
  Point<3, int> point4{-1, -2, -3};
  Point<3, int> point5{3, 6, 9};
};

/*
 * Demonstrate some basic vector arithmetic
 */
TEST_F(BezierTestingSuite, PointTests) {
  // Expect equality.
  EXPECT_EQ(point1, point2 + point3);
  EXPECT_EQ(point1 - point2, point3);
  EXPECT_EQ((-1) * point1, point4);
  EXPECT_EQ(point1 * 3, point5);
  EXPECT_EQ(-point1, point4);
}

TEST_F(BezierTestingSuite, PointCreationTests) {
  EXPECT_TRUE((utils::type_traits::isPoint_v<Point<3, double>>));
  EXPECT_FALSE((utils::type_traits::isPoint_v<BezierSpline<2,Point<3, double>, double>>));
}

}  // namespace bezman::tests::point_test