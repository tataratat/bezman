#include <cmath>

#include "bezierManipulation/src/bezier_spline_group.hpp"
#include "bezierManipulation/src/utils/export.hpp"

using namespace beziermanipulation;

// Define aliases
using Point2D = Point<2, double>;
using BezierGroup = BezierSplineGroup<2, Point2D, double>;
using Bezier = BezierSpline<2, Point2D, double>;

BezierGroup SimpleCrossTile() {
  // Degrees
  std::array<std::size_t, 2> degrees{1, 1};
  std::array<std::size_t, 2> second_order_degrees{1, 2};

  // Auxiliary values
  const double one_third = 1. / 3.;

  // CTPS
  std::vector<Point2D> ctps_center{
      Point2D{one_third, 0.25}, Point2D{2. * one_third, 0.25},
      Point2D{one_third, 0.5}, Point2D{2. * one_third, 0.5}};

  std::vector<Point2D> ctps_left{
      Point2D{one_third, 0.25}, Point2D{one_third, 0.5}, Point2D{0., 0.25},
      Point2D{0., 0.5},         Point2D{0., 0.0},       Point2D{0., 1.}};

  std::vector<Point2D> ctps_right{Point2D{2. * one_third, 0.5},
                                  Point2D{2. * one_third, 0.25},
                                  Point2D{1., 0.5},
                                  Point2D{1., 0.25},
                                  Point2D{1., 1.},
                                  Point2D{1., 0.0}};

  std::vector<Point2D> ctps_down{Point2D{2. * one_third, 0.25},
                                 Point2D{one_third, 0.25}, Point2D{1.0, 0.0},
                                 Point2D{0.0, 0.0}};

  std::vector<Point2D> ctps_up{
      Point2D{one_third, 0.5}, Point2D{2. * one_third, 0.5},
      Point2D{0., 0.85},       Point2D{1., 0.85},
      Point2D{0., 1.},         Point2D{1., 1.}};

  return BezierGroup{
      Bezier{degrees, ctps_center}, Bezier{second_order_degrees, ctps_left},
      Bezier{second_order_degrees, ctps_right}, Bezier{degrees, ctps_down},
      Bezier{second_order_degrees, ctps_up}};
}
BezierGroup SimpleCrossTile(const double thickness, const bool deriv = false) {
  std::array<std::size_t, 2> degrees{1, 1};
  if (!deriv) {
    std::vector<Point2D> ctps_center{
        Point2D{0.5 - 0.5 * thickness, 0.5 - 0.5 * thickness},
        Point2D{0.5 + 0.5 * thickness, 0.5 - 0.5 * thickness},
        Point2D{0.5 - 0.5 * thickness, 0.5 + 0.5 * thickness},
        Point2D{0.5 + 0.5 * thickness, 0.5 + 0.5 * thickness}};
    std::vector<Point2D> ctps_left{
        Point2D{0., 0.5 - 0.5 * thickness},
        Point2D{0.5 - 0.5 * thickness, 0.5 - 0.5 * thickness},
        Point2D{0., 0.5 + 0.5 * thickness},
        Point2D{0.5 - 0.5 * thickness, 0.5 + 0.5 * thickness}};
    std::vector<Point2D> ctps_right{
        Point2D{0.5 + 0.5 * thickness, 0.5 - 0.5 * thickness},
        Point2D{1., 0.5 - 0.5 * thickness},
        Point2D{0.5 + 0.5 * thickness, 0.5 + 0.5 * thickness},
        Point2D{1., 0.5 + 0.5 * thickness}};
    std::vector<Point2D> ctps_down{
        Point2D{0.5 - 0.5 * thickness, 0.}, Point2D{0.5 + 0.5 * thickness, 0.},
        Point2D{0.5 - 0.5 * thickness, 0.5 - 0.5 * thickness},
        Point2D{0.5 + 0.5 * thickness, 0.5 - 0.5 * thickness}};
    std::vector<Point2D> ctps_up{
        Point2D{0.5 - 0.5 * thickness, 0.5 + 0.5 * thickness},
        Point2D{0.5 + 0.5 * thickness, 0.5 + 0.5 * thickness},
        Point2D{0.5 - 0.5 * thickness, 1.}, Point2D{0.5 + 0.5 * thickness, 1.}};

    return BezierGroup{Bezier{degrees, ctps_center}, Bezier{degrees, ctps_left},
                       Bezier{degrees, ctps_right}, Bezier{degrees, ctps_down},
                       Bezier{degrees, ctps_up}};
  } else {
    std::vector<Point2D> ctps_center_deriv{
        Point2D{-0.5, -0.5}, Point2D{0.5, -0.5}, Point2D{-0.5, 0.5},
        Point2D{0.5, 0.5}};
    std::vector<Point2D> ctps_left_deriv{Point2D{0., -0.5}, Point2D{-0.5, -0.5},
                                         Point2D{0., 0.5}, Point2D{-0.5, 0.5}};
    std::vector<Point2D> ctps_right_deriv{Point2D{0.5, -0.5}, Point2D{0., -0.5},
                                          Point2D{0.5, 0.5}, Point2D{0., 0.5}};
    std::vector<Point2D> ctps_down_deriv{Point2D{-0.5, 0.}, Point2D{0.5, 0.},
                                         Point2D{-0.5, -0.5},
                                         Point2D{0.5, -0.5}};
    std::vector<Point2D> ctps_up_deriv{Point2D{-0.5, 0.5}, Point2D{0.5, 0.5},
                                       Point2D{-0.5, 0.}, Point2D{0.5, 0.}};

    return BezierGroup{
        Bezier{degrees, ctps_center_deriv}, Bezier{degrees, ctps_left_deriv},
        Bezier{degrees, ctps_right_deriv}, Bezier{degrees, ctps_down_deriv},
        Bezier{degrees, ctps_up_deriv}};
  }
}

BezierGroup BulkSquare() {
  // Create Crosstile Group
  std::vector<Point2D> ctps_deformation_function_bulk{
      Point2D{1., 0.},   Point2D{1.5, -0.2}, Point2D{2., 0.},
      Point2D{0.8, 0.5}, Point2D{1.5, 0.5},  Point2D{2.2, 0.5},
      Point2D{1., 1.},   Point2D{1.5, 1.2},  Point2D{2., 1.}};
  std::array<std::size_t, 2> deformation_function_degrees{2, 2};
  return BezierGroup{
      Bezier{deformation_function_degrees, ctps_deformation_function_bulk}};
}

BezierGroup CircleGroup(const double innerR, const double outerR,
                        const int segments, const double arc_degrees) {
  ///
  BezierGroup ringsegments{segments};
  constexpr double PI = acos(-1);
  const double degrees_per_segment = arc_degrees / segments;
  const double outdor_fact = 1. / std::sin(PI / 2. - degrees_per_segment / 2.);
  const std::array<std::size_t, 2> degrees{2, 1};
  for (int i_segment{}; i_segment < segments; i_segment++) {
    const double startsin = std::sin(degrees_per_segment * i_segment);
    const double startcos = std::cos(degrees_per_segment * i_segment);
    const double middlesin = std::sin(degrees_per_segment * (i_segment + .5));
    const double middlecos = std::cos(degrees_per_segment * (i_segment + .5));
    const double endsin = std::sin(degrees_per_segment * (i_segment + 1.));
    const double endcos = std::cos(degrees_per_segment * (i_segment + 1.));
    std::vector<Point2D> ctps{Point2D{startcos * outerR, startsin * outerR},
                              Point2D{middlecos * outdor_fact * outerR,
                                      middlesin * outdor_fact * outerR},
                              Point2D{endcos * outerR, endsin * outerR},
                              Point2D{startcos * innerR, startsin * innerR},
                              Point2D{middlecos * outdor_fact * innerR,
                                      middlesin * outdor_fact * innerR},
                              Point2D{endcos * innerR, endsin * innerR}};

    ringsegments[i_segment] =
        Bezier{degrees, ctps}.OrderElevateAlongParametricDimension(1);
  }
  return ringsegments;
}

int main() {
  const int n_segments = 20;
  const double thickness = 0.3;
  // Inner function
  auto microtile = SimpleCrossTile(thickness);

  // Outer Function
  auto deformation_function =
      CircleGroup(1., 2., n_segments, 2*std::acos(-1));

  // Compose composition
  const auto test_composition = deformation_function.Compose(microtile);

  utils::Export::AsJSON(test_composition, "composed_microstructure");
  utils::Export::AsMFEM(test_composition, "composed_microstructure");
  return 0;
}
