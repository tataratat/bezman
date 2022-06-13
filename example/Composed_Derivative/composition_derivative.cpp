// Load Bezier
#include "bezman/src/bezier_spline_group.hpp"
#include "bezman/src/utils/export.hpp"

using namespace bezman;

int main() {
  // Define aliases
  using Point2D = Point<2, double>;

  // Create 2D CrossTile CTPS
  std::vector<Point2D> ctps_center{Point2D{0.2, 0.}, Point2D{0.4, 0.},
                                   Point2D{0.0, 1.}, Point2D{1.0, 1.}};
  std::vector<Point2D> ctps_deformation_function_bulk{
      Point2D{1., 0.},   Point2D{1.5, -0.2}, Point2D{2., 0.},
      Point2D{0.8, 0.5}, Point2D{1.5, 0.5},  Point2D{2.2, 0.5},
      Point2D{2., 2.},   Point2D{2.5, 1.8},  Point2D{3., 1.}};

  std::array<std::size_t, 2> cross_tile_degrees{1, 1};
  std::array<std::size_t, 2> deformation_function_degrees{2, 2};

  // Create Crosstile Group
  BezierSpline<2, Point2D, double> microtile_cross{cross_tile_degrees,
                                                   ctps_center};

  // Create outer Spline
  BezierSpline<2, Point2D, double> deformation_function =
      BezierSpline<2, Point2D, double>(deformation_function_degrees,
                                       ctps_deformation_function_bulk);

  // Derivative of the inner function
  BezierSpline<2, Point2D, double> microtile_cross_derivative{microtile_cross};

  // Set all points to zero
  for (std::size_t i_point{};
       i_point < microtile_cross_derivative.control_points.size(); i_point++) {
    microtile_cross_derivative.control_points[i_point] = Point2D{0., 0.};
  }
  // Define the derivative of the inner function (Here just moving the upper
  // most points to the outside)
  microtile_cross_derivative.control_points[2] = Point2D{-1., 0.};
  microtile_cross_derivative.control_points[3] = Point2D{+1., 0.};

  auto microstructure_derivative =
      deformation_function.DerivativeWRTParametricDimension(0).Compose(
          microtile_cross) *
      microtile_cross_derivative.ExtractDimension(0);
  microstructure_derivative +=
      deformation_function.DerivativeWRTParametricDimension(1).Compose(
          microtile_cross) *
      microtile_cross_derivative.ExtractDimension(1);

  const auto test_composition = deformation_function.Compose(microtile_cross);
  utils::Export::GuessByExtension(test_composition,
                                  "composed_microstructure.xml");
  utils::Export::GuessByExtension(microstructure_derivative,
                                  "composed_microstructure_derivative.xml");
  utils::Export::GuessByExtension(microtile_cross_derivative,
                                  "inner_function_derived.xml");
  utils::Export::GuessByExtension(microtile_cross, "inner_function.xml");
  utils::Export::GuessByExtension(deformation_function, "outer_function.xml");
  return 0;
}
