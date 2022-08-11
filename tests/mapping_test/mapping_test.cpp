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

namespace bezman::tests::mapping_test {

class BezierTestingSuite : public ::testing::Test {
  using Point3D = Point<3, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> surface_ctps{Point3D{-1., 0., 2.5}, Point3D{2., 1., 2.5},
                                    Point3D{-1., 0., 4.0},
                                    Point3D{2., 1., 4.0}};
  std::vector<Point3D> ref_surface_ctps{
      Point3D{0., 0., 0.}, Point3D{1., 1., 0.}, Point3D{0., 0., 1.},
      Point3D{1., 1., 1.}};

  std::array<std::size_t, 2> surface_degrees{1, 1};

  // Define Splines
  BezierSpline<2, Point3D, double> surface =
      BezierSpline<2, Point3D, double>(surface_degrees, surface_ctps);
  BezierSpline<2, Point3D, double> reference_surface =
      BezierSpline<2, Point3D, double>(surface_degrees, ref_surface_ctps);
};

/*
 * Check the transposition functions that are all combined in the fit to unit
 * cube option
 */
TEST_F(BezierTestingSuite, MappingToUnitCube) {
  // Expect equality.
  EXPECT_FALSE(surface.FitsIntoUnitCube());
  surface.FitIntoUnitCube();
  EXPECT_TRUE(surface.FitsIntoUnitCube());
  EXPECT_TRUE(reference_surface.FitsIntoUnitCube());
  for (std::size_t i{}; i < surface.GetNumberOfControlPoints(); i++) {
    for (std::size_t i_dim{}; i_dim < 3; i_dim++) {
      EXPECT_FLOAT_EQ(surface.control_points[i][i_dim],
                      ref_surface_ctps[i][i_dim]);
    }
  }
}

}  // namespace bezman::tests::mapping_test
