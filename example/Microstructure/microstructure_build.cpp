#include <cmath>

#include "bezman/src/bezier_spline_group.hpp"
#include "bezman/src/utils/export.hpp"

using namespace bezman;

// Define aliases
using Point2D = Point<2, double>;
using BezierGroup = BezierSplineGroup<2, Point2D, double>;
using Bezier = BezierSpline<2, Point2D, double>;

/**
 * @brief Creates a simple linear crosstile with constant thickness
 */
BezierGroup SimpleCrossTile(const double thickness) {
  std::array<std::size_t, 2> degrees{1, 1};
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
}

/**
 * @brief Creates a segmented ring to be filled with tiles
 */
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

    ringsegments[i_segment] = Bezier{degrees, ctps};
  }
  return ringsegments;
}

int main() {
  // Parameters for the outer function
  const int n_segments = 20;
  const double innerRadius{1.}, outerRadius{2.};
  const double arc_segment{2 * std::acos(-1)};

  // Microtile Thickness
  const double thickness = 0.3;

  // Inner function
  auto microtile = SimpleCrossTile(thickness);

  // Outer Function
  auto deformation_function =
      CircleGroup(innerRadius, outerRadius, n_segments, arc_segment);

  // Compose composition
  const auto test_composition = deformation_function.Compose(microtile);

  // Export the composed structure in different Formats for testing
  utils::Export::AsJSON(test_composition, "composed_microstructure");
  utils::Export::AsMFEM(test_composition, "composed_microstructure");
  utils::Export::AsIRIT(test_composition, "composed_microstructure");
  utils::Export::AsXML(test_composition, "composed_microstructure");
  return 0;
}
