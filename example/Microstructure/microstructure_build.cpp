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

#include <cmath>

#include "bezman/src/bezier_group.hpp"
#include "bezman/src/rational_bezier_spline.hpp"
#include "bezman/src/utils/export.hpp"

using namespace bezman;

// Define aliases
using Point2D = Point<2, double>;
using PolynomialBezier = BezierSpline<2, Point2D, double>;
using RationalBezier = RationalBezierSpline<2, Point2D, double>;
using PolynomialBezierGroup = BezierGroup<PolynomialBezier>;
using RationalBezierGroup = BezierGroup<RationalBezier>;

/**
 * @brief Creates a simple linear crosstile with constant thickness
 */
PolynomialBezierGroup SimpleCrossTile(const double thickness) {
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

  return PolynomialBezierGroup{PolynomialBezier{degrees, ctps_center},
                               PolynomialBezier{degrees, ctps_left},
                               PolynomialBezier{degrees, ctps_right},
                               PolynomialBezier{degrees, ctps_down},
                               PolynomialBezier{degrees, ctps_up}};
}

/**
 * @brief Creates a segmented ring to be filled with tiles
 */
RationalBezierGroup CircleGroup(const double innerR, const double outerR,
                                const int segments, const double arc_degrees) {
  ///
  RationalBezierGroup ringsegments{segments};
  constexpr double PI = acos(-1);
  const double degrees_per_segment = arc_degrees / segments;
  const double single_weight{std::sin(PI / 2. - degrees_per_segment / 2.)};
  const double outdor_fact = 1. / single_weight;
  const std::array<std::size_t, 2> degrees{2, 1};
  std::vector<double> weights{1, single_weight, 1, 1, single_weight, 1};
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

    ringsegments[i_segment] = RationalBezier{degrees, ctps, weights};
  }
  return ringsegments;
}

int main() {
  // Parameters for the outer function
  const int n_segments = 4;
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
  return 0;
}
