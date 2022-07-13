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

  // Points for testing division
  Point2D point_div{2., 5.};
  Point2D point_div_inv_x2{1., 0.4};
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

  EXPECT_FLOAT_EQ((point_div / 2.)[0], 1.);
  EXPECT_FLOAT_EQ((point_div / 2.)[1], 2.5);
  EXPECT_FLOAT_EQ((2. / point_div)[0], point_div_inv_x2[0]);
  EXPECT_FLOAT_EQ((2. / point_div)[1], point_div_inv_x2[1]);
  point_div /= 5.;
  EXPECT_FLOAT_EQ(point_div[0], 0.4);
  EXPECT_FLOAT_EQ(point_div[1], 1.);
}

TEST_F(BezierTestingSuite, PointCreationTests) {
  EXPECT_TRUE((utils::type_traits::isPoint_v<Point<3, double>>));
  EXPECT_FALSE((utils::type_traits::isPoint_v<
                BezierSpline<2, Point<3, double>, double>>));
}

}  // namespace bezman::tests::point_test