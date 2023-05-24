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
#include "bezman/src/rational_bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::composition_test {

class BezierTestingSuite : public ::testing::Test {
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // void TearDown() override {}

  // Provide some data to form splines.
  std::vector<Point3D> surface_ctps{Point3D{0., 0., 0.}, Point3D{2., 1., 0.},
                                    Point3D{0., 2., 0.}, Point3D{2., 3., 0.}};
  std::vector<Point3D> surface_ctps2{Point3D{0., 0., 0.}, Point3D{1., 1., 0.},
                                     Point3D{2., 1., 0.}, Point3D{0., 2., 0.},
                                     Point3D{1., 2., 0.}, Point3D{2., 3., 0.}};
  std::vector<Point2D> line_ctps{Point2D{0., 0.}, Point2D{1., 0.},
                                 Point2D{1., 1.}};
  std::vector<double> line_weights{1., 1 / std::sqrt(2), 1.};
  std::vector<double> surface_weights{1., 1 / std::sqrt(2), 1.,
                                      1., 1 / std::sqrt(2), 1.};

  std::array<std::size_t, 2> surface_degrees{1, 1};
  std::array<std::size_t, 2> surface_degrees2{2, 1};
  std::array<std::size_t, 1> line_degrees{2};

  // Define Splines
  BezierSpline<2, Point3D, double> surface =
      BezierSpline<2, Point3D, double>(surface_degrees, surface_ctps);
  RationalBezierSpline<2, Point3D, double> r_surface =
      RationalBezierSpline<2, Point3D, double>(surface_degrees2, surface_ctps2,
                                               surface_weights);
  BezierSpline<1, Point2D, double> line =
      BezierSpline<1, Point2D, double>(line_degrees, line_ctps);
  RationalBezierSpline<1, Point2D, double> r_line =
      RationalBezierSpline<1, Point2D, double>(line_degrees, line_ctps,
                                               line_weights);
};

/*
 * Demonstrate spline composition sensitivities Polynomial(Polynomial)
 *
 * This is done using finite differences for validation and comparing the
 * approximate solution to the analytical one.
 */
TEST_F(BezierTestingSuite, CompositionSensitivityPolyPoly) {
  // Conpute baseline values
  const auto composed_spline = surface.Compose(line);
  const auto composed_sensitivity = surface.ComposeSensitivity(line);

  // Baseline check
  EXPECT_EQ(surface.GetNumberOfControlPoints(), composed_sensitivity.size());
  for (std::size_t i{}; i < surface.GetNumberOfControlPoints(); i++) {
    EXPECT_EQ(composed_spline.GetNumberOfControlPoints(),
              composed_sensitivity[i].GetNumberOfControlPoints());
  }

  // Set finite difference step width
  const double dx{1e-5}, inv_dx{1. / dx};
  const double tolerance{1e-5};

  // Test the control point positions against the finite element solution for
  // every control point in the outer spline (in all physical directions)
  for (std::size_t i{}; i < surface.GetNumberOfControlPoints(); i++) {
    for (std::size_t j{}; j < 2; j++) {
      // Move control point
      surface.control_points[i][j] += dx;

      // Perform composition
      const auto composed_spline_dx = surface.Compose(line);

      // Check values
      for (std::size_t k{}; k < surface.GetNumberOfControlPoints(); k++) {
        for (std::size_t l{}; l < 2; l++) {
          if (l == j) {
            const double fd_value = (composed_spline_dx.control_points[k][l] -
                                     composed_spline.control_points[k][l]) *
                                    inv_dx;

            // Differentiation is a linear operation in this case, so we can
            // expect the values to be equal
            EXPECT_TRUE(
                std::abs(fd_value - composed_sensitivity[i].control_points[k]) <
                tolerance);
          } else {
            EXPECT_DOUBLE_EQ(composed_spline_dx.control_points[k][l],
                             composed_spline.control_points[k][l]);
          }
        }
      }
      // Clean up - move control point back
      surface.control_points[i][j] -= dx;
    }
  }
}

/*
 * Demonstrate spline composition sensitivities Polynomial(Rational)
 *
 * This is done using finite differences for validation and comparing the
 * approximate solution to the analytical one.
 */
TEST_F(BezierTestingSuite, CompositionSensitivityPolyRat) {
  // Conpute baseline values
  const auto composed_spline = surface.Compose(r_line);
  const auto composed_sensitivity = surface.ComposeSensitivity(r_line);

  // Baseline check
  EXPECT_EQ(surface.GetNumberOfControlPoints(), composed_sensitivity.size());
  for (std::size_t i{}; i < surface.GetNumberOfControlPoints(); i++) {
    EXPECT_EQ(composed_spline.GetNumberOfControlPoints(),
              composed_sensitivity[i].GetNumberOfControlPoints());
  }

  // Set finite difference step width
  const double dx{1e-5}, inv_dx{1. / dx};
  const double tolerance{1e-5};

  // Test the control point positions against the finite element solution for
  // every control point in the outer spline (in all physical directions)
  for (std::size_t i{}; i < surface.GetNumberOfControlPoints(); i++) {
    for (std::size_t j{}; j < 2; j++) {
      // Move control point
      surface.control_points[i][j] += dx;

      // Perform composition
      const auto composed_spline_dx = surface.Compose(r_line);

      // Check values
      for (std::size_t k{}; k < surface.GetNumberOfControlPoints(); k++) {
        for (std::size_t l{}; l < 2; l++) {
          if (l == j) {
            const double fd_value =
                (composed_spline_dx.GetWeightedControlPoints()[k][l] -
                 composed_spline.GetWeightedControlPoints()[k][l]) *
                inv_dx;
            EXPECT_TRUE(
                std::abs(
                    fd_value -
                    composed_sensitivity[i].GetWeightedControlPoints()[k]) <
                tolerance);
          } else {
            EXPECT_DOUBLE_EQ(
                composed_spline_dx.GetWeightedControlPoints()[k][l],
                composed_spline.GetWeightedControlPoints()[k][l]);
          }
        }
        EXPECT_DOUBLE_EQ(composed_spline_dx.GetWeights()[k],
                         composed_spline.GetWeights()[k]);
      }
      // Clean up - move control point back
      surface.control_points[i][j] -= dx;
    }
  }
}

/*
 * Demonstrate spline composition sensitivities Rational(Polynomial)
 *
 * This is done using finite differences for validation and comparing the
 * approximate solution to the analytical one.
 */
TEST_F(BezierTestingSuite, CompositionSensitivityRatPoly) {
  // Conpute baseline values
  const auto composed_spline = r_surface.Compose(line);
  const auto composed_sensitivity = r_surface.ComposeSensitivity(line);

  // Baseline check
  EXPECT_EQ(r_surface.GetNumberOfControlPoints(), composed_sensitivity.size());
  for (std::size_t i{}; i < r_surface.GetNumberOfControlPoints(); i++) {
    EXPECT_EQ(composed_spline.GetNumberOfControlPoints(),
              composed_sensitivity[i].GetNumberOfControlPoints());
  }

  // Set finite difference step width
  const double dx{1e-5}, inv_dx{1. / dx};
  const double tolerance{1e-5};

  // Test the control point positions against the finite element solution for
  // every control point in the outer spline (in all physical directions)
  for (std::size_t i{}; i < r_surface.GetNumberOfControlPoints(); i++) {
    const double dx_local = dx * r_surface.GetWeights()[i];
    for (std::size_t j{}; j < 2; j++) {
      // Move control point
      r_surface.GetWeightedControlPoints()[i][j] += dx_local;

      // Perform composition
      const auto composed_spline_dx = r_surface.Compose(line);

      // Check values
      for (std::size_t k{}; k < r_surface.GetNumberOfControlPoints(); k++) {
        for (std::size_t l{}; l < 2; l++) {
          if (l == j) {
            const double fd_value =
                (composed_spline_dx.GetWeightedControlPoints()[k][l] -
                 composed_spline.GetWeightedControlPoints()[k][l]) *
                inv_dx;
            EXPECT_TRUE(
                std::abs(
                    fd_value -
                    composed_sensitivity[i].GetWeightedControlPoints()[k]) <
                tolerance);
          } else {
            EXPECT_DOUBLE_EQ(
                composed_spline_dx.GetWeightedControlPoints()[k][l],
                composed_spline.GetWeightedControlPoints()[k][l]);
          }
        }
        EXPECT_DOUBLE_EQ(composed_spline_dx.GetWeights()[k],
                         composed_spline.GetWeights()[k]);
      }
      // Clean up - move control point back
      r_surface.GetWeightedControlPoints()[i][j] -= dx_local;
    }
  }
}

/*
 * Demonstrate spline composition sensitivities Rational(Rational)
 *
 * This is done using finite differences for validation and comparing the
 * approximate solution to the analytical one.
 */
TEST_F(BezierTestingSuite, CompositionSensitivityRatRat) {
  // Conpute baseline values
  const auto composed_spline = r_surface.Compose(r_line);
  const auto composed_sensitivity = r_surface.ComposeSensitivity(r_line);

  // Baseline check
  EXPECT_EQ(r_surface.GetNumberOfControlPoints(), composed_sensitivity.size());
  for (std::size_t i{}; i < r_surface.GetNumberOfControlPoints(); i++) {
    EXPECT_EQ(composed_spline.GetNumberOfControlPoints(),
              composed_sensitivity[i].GetNumberOfControlPoints());
  }

  // Set finite difference step width
  const double dx{1e-5}, inv_dx{1. / dx};
  const double tolerance{1e-5};

  // Test the control point positions against the finite element solution for
  // every control point in the outer spline (in all physical directions)
  for (std::size_t i{}; i < r_surface.GetNumberOfControlPoints(); i++) {
    const double dx_local = dx * r_surface.GetWeights()[i];
    for (std::size_t j{}; j < 2; j++) {
      // Move control point
      r_surface.GetWeightedControlPoints()[i][j] += dx_local;

      // Perform composition
      const auto composed_spline_dx = r_surface.Compose(r_line);

      // Check values
      for (std::size_t k{}; k < r_surface.GetNumberOfControlPoints(); k++) {
        for (std::size_t l{}; l < 2; l++) {
          if (l == j) {
            const double fd_value =
                (composed_spline_dx.GetWeightedControlPoints()[k][l] -
                 composed_spline.GetWeightedControlPoints()[k][l]) *
                inv_dx;
            EXPECT_TRUE(
                std::abs(
                    fd_value -
                    composed_sensitivity[i].GetWeightedControlPoints()[k]) <
                tolerance);
          } else {
            EXPECT_DOUBLE_EQ(
                composed_spline_dx.GetWeightedControlPoints()[k][l],
                composed_spline.GetWeightedControlPoints()[k][l]);
          }
        }
        EXPECT_DOUBLE_EQ(composed_spline_dx.GetWeights()[k],
                         composed_spline.GetWeights()[k]);
      }
      // Clean up - move control point back
      r_surface.GetWeightedControlPoints()[i][j] -= dx_local;
    }
  }
}

}  // namespace bezman::tests::composition_test