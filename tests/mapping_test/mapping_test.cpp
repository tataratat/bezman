#include <gtest/gtest.h>

#include <array>

#define MAX_BINOMIAL_DEGREE 62u  // Required for second test set
#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

namespace beziermanipulation::tests::mapping_test {

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
  surface.FitToUnitCube();
  EXPECT_TRUE(surface.FitsIntoUnitCube());
  EXPECT_TRUE(reference_surface.FitsIntoUnitCube());
  for (std::size_t i{}; i < surface.NumberOfControlPoints; i++) {
    for (std::size_t i_dim{}; i_dim < 3; i_dim++) {
      EXPECT_FLOAT_EQ(surface.control_points[i][i_dim],
                      ref_surface_ctps[i][i_dim]);
    }
  }
}

}  // namespace beziermanipulation::tests::mapping_test