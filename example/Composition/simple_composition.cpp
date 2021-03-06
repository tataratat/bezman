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

// Load Bezier
#include "bezman/src/bezier_spline_group.hpp"
#include "bezman/src/utils/export.hpp"

using namespace bezman;

int main() {
  // Define aliases
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

  // Create 2D CrossTile CTPS
  const double one_third{1. / 3.};
  std::vector<Point2D> ctps_center{
      Point2D{one_third, one_third}, Point2D{2 * one_third, one_third},
      Point2D{one_third, 2 * one_third}, Point2D{2 * one_third, 2 * one_third}};
  std::vector<Point2D> ctps_deformation_function_slim{
      Point2D{0., 0.},   Point2D{0.5, 0.2}, Point2D{1., 0.},
      Point2D{0.2, 0.5}, Point2D{0.5, 0.5}, Point2D{.8, 0.5},
      Point2D{0., 1.},   Point2D{0.5, 0.8}, Point2D{1., 1.}};
  std::vector<Point2D> ctps_deformation_function_bulk{
      Point2D{1., 0.},   Point2D{1.5, -0.2}, Point2D{2., 0.},
      Point2D{0.8, 0.5}, Point2D{1.5, 0.5},  Point2D{2.2, 0.5},
      Point2D{1., 1.},   Point2D{1.5, 1.2},  Point2D{2., 1.}};

  std::array<std::size_t, 2> cross_tile_degrees{1, 1};
  std::array<std::size_t, 2> deformation_function_degrees{2, 2};

  // Create Crosstile Group
  BezierSpline<2, Point2D, double> center_surface{cross_tile_degrees,
                                                  ctps_center};
  BezierSplineGroup<2, Point2D, double> microtile_cross{
      center_surface, center_surface + Point2D{one_third, 0.},
      center_surface + Point2D{0., one_third},
      center_surface - Point2D{one_third, 0.},
      center_surface - Point2D{0., one_third}};
  // Create outer Spline
  BezierSplineGroup<2, Point2D, double> deformation_function{
      BezierSpline<2, Point2D, double>(deformation_function_degrees,
                                       ctps_deformation_function_slim),
      BezierSpline<2, Point2D, double>(deformation_function_degrees,
                                       ctps_deformation_function_bulk)};

  const auto test_composition = deformation_function.Compose(microtile_cross);
  utils::Export::GuessByExtension(test_composition,
                                  "composed_microstructure.xml");
  return 0;
}
