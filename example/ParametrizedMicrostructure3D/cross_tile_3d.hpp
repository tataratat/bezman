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

#ifndef EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROTILE_HPP
#define EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROTILE_HPP

#include <algorithm>
#include <vector>

#include "bezman/src/bezier_group.hpp"
#include "bezman/src/utils/computational_differentiation/algo_diff_type.hpp"

using namespace bezman;

/**
 * @brief Cross tile with variable thickness
 *
 * Defines a Microtile-generator, that creates tiles based on the evaluation of
 * a function at specific evaluation points (defined on the unit cube). These
 * points are set at the edges of the tile's branches
 *
 * The branches are defined in a way that their attachement points are vertical
 * to the unit cube surface.
 * This class works with doubles as basis type.
 */
class CrossTile3D {
 private:
  // Aliases
  using ADT = utils::computational_differentiation::AlgoDiffType<double>;
  using PointADT3D = Point<3, ADT>;
  using Point3D = Point<3, double>;
  using Bezier = BezierSpline<3, Point3D, double>;
  using PolyBezierGroup = BezierGroup<Bezier>;
  /**
   * @brief Extracts the Derivatives from the derivative vector of each control
   * point and defines new sets of splines, which  are then used to determine
   * the derivatives of the microtiles
   */
  static std::vector<Point3D> ExtractDerivativeFromPointCloud(
      const std::vector<PointADT3D>& points, const std::size_t& derivative) {
    // Initialize return type
    std::vector<Point3D> return_vector(points.size());

    for (std::size_t i_point{}; i_point < points.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < PointADT3D::kSpatialDimension;
           i_dim++) {
        return_vector[i_point][i_dim] =
            points[i_point][i_dim].GetDerivatives()[derivative];
      }
    }
    return return_vector;
  }

  /**
   * @brief Extracts the values from each control point to change the type of
   * the points.
   *
   * @param points
   * @return std::vector<Point3D>
   */
  static std::vector<Point3D> ExtractValuesFromPointCloud(
      const std::vector<PointADT3D>& points) {
    // Initialize return type
    std::vector<Point3D> return_vector(points.size());
    for (std::size_t i_point{}; i_point < points.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < PointADT3D::kSpatialDimension;
           i_dim++) {
        return_vector[i_point][i_dim] = points[i_point][i_dim].GetValue();
      }
    }
    return return_vector;
  }

 public:
  /// Number of Splines per Tile
  static constexpr unsigned int kNumberOfSplines = 7;
  /// Number of Evaluation Points
  static constexpr std::size_t kNumberOfEvaluationPoints = 7;
  /// Evaluation Points
  static constexpr std::array<Point3D, kNumberOfEvaluationPoints>
      kEvaluationPoints{Point3D{.5, 0.5, 0.5}, Point3D{.5, 0., 0.5},
                        Point3D{1., 0.5, 0.5}, Point3D{.5, 1., 0.5},
                        Point3D{0., 0.5, 0.5}, Point3D{0.5, 0.5, 0.},
                        Point3D{0.5, 0.5, 1.}};

  /// Generator Function for Derivatives
  static std::vector<PolyBezierGroup> GenerateMicrostructureDerivatives(
      const std::array<ADT, kNumberOfEvaluationPoints>& evaluations) {
    // Check if all evaluations have the same number of derivatives
    const auto number_of_derivatives = evaluations[0].GetNumberOfDerivatives();
    assert(std::all_of(evaluations.begin(), evaluations.end(),
                       [&number_of_derivatives](const ADT& ev_point) -> bool {
                         return (ev_point.GetNumberOfDerivatives() ==
                                 number_of_derivatives);
                       }));
    // Precomputed values and aliases
    const double one_third{1. / 3.};
    const double two_thirds{2. / 3.};

    // Name abreviations (only readability no purpose)
    const auto& center_width = evaluations[0];
    const auto& thickness_low = evaluations[1];
    const auto& thickness_right = evaluations[2];
    const auto& thickness_up = evaluations[3];
    const auto& thickness_left = evaluations[4];
    const auto& thickness_front = evaluations[5];
    const auto& thickness_back = evaluations[6];

    // Define degrees
    std::array<std::size_t, 3> center_degrees{1, 1, 1};
    std::array<std::size_t, 3> horizontal_degrees{3, 1, 1};
    std::array<std::size_t, 3> vertical_degrees{1, 3, 1};
    std::array<std::size_t, 3> pass_z_degrees{1, 1, 3};

    /*
     *              *--------*
     *              |   up   |
     *       *------*--------*-------*
     *       | left | center | right |
     *       *------*--------*-------*
     *              |  down  |
     *              *--------*
     *
     * Definition of the control points.
     * order is required to get 90 deg angles from center and from start of
     * individual branches (center point tangents might not be required)
     *
     */

    // center ctps
    std::vector<PointADT3D> ctps_center{
        PointADT3D{0.5 * (1. - center_width), 0.5 * (1. - center_width),
                   0.5 * (1. - center_width)},
        PointADT3D{0.5 * (1. + center_width), 0.5 * (1. - center_width),
                   0.5 * (1. - center_width)},
        PointADT3D{0.5 * (1. - center_width), 0.5 * (1. + center_width),
                   0.5 * (1. - center_width)},
        PointADT3D{0.5 * (1. + center_width), 0.5 * (1. + center_width),
                   0.5 * (1. - center_width)},
        PointADT3D{0.5 * (1. - center_width), 0.5 * (1. - center_width),
                   0.5 * (1. + center_width)},
        PointADT3D{0.5 * (1. + center_width), 0.5 * (1. - center_width),
                   0.5 * (1. + center_width)},
        PointADT3D{0.5 * (1. - center_width), 0.5 * (1. + center_width),
                   0.5 * (1. + center_width)},
        PointADT3D{0.5 * (1. + center_width), 0.5 * (1. + center_width),
                   0.5 * (1. + center_width)}};

    // Lower CTPS
    std::vector<PointADT3D> ctps_low{
        PointADT3D{0.5 * (1. - thickness_low), ADT(0., number_of_derivatives),
                   0.5 * (1. - thickness_low)},
        PointADT3D{0.5 * (1. + thickness_low), ADT(0., number_of_derivatives),
                   0.5 * (1. - thickness_low)},
        PointADT3D{0.5 * (1. - thickness_low), one_third * ctps_center[0][1],
                   0.5 * (1. - thickness_low)},
        PointADT3D{0.5 * (1. + thickness_low), one_third * ctps_center[0][1],
                   0.5 * (1. - thickness_low)},
        PointADT3D{ctps_center[0][0], two_thirds * ctps_center[0][1],
                   ctps_center[0][2]},
        PointADT3D{ctps_center[1][0], two_thirds * ctps_center[0][1],
                   ctps_center[1][2]},
        ctps_center[0],
        ctps_center[1],
        PointADT3D{0.5 * (1. - thickness_low), ADT(0., number_of_derivatives),
                   0.5 * (1. + thickness_low)},
        PointADT3D{0.5 * (1. + thickness_low), ADT(0., number_of_derivatives),
                   0.5 * (1. + thickness_low)},
        PointADT3D{0.5 * (1. - thickness_low), one_third * ctps_center[0][1],
                   0.5 * (1. + thickness_low)},
        PointADT3D{0.5 * (1. + thickness_low), one_third * ctps_center[0][1],
                   0.5 * (1. + thickness_low)},
        PointADT3D{ctps_center[4][0], two_thirds * ctps_center[0][1],
                   ctps_center[4][2]},
        PointADT3D{ctps_center[5][0], two_thirds * ctps_center[0][1],
                   ctps_center[5][2]},
        ctps_center[4],
        ctps_center[5]};

    // Upper CTPS
    std::vector<PointADT3D> ctps_up{
        ctps_center[2],
        ctps_center[3],
        PointADT3D{ctps_center[2][0], one_third * (2. * ctps_center[2][1] + 1.),
                   ctps_center[2][2]},
        PointADT3D{ctps_center[3][0], one_third * (2. * ctps_center[2][1] + 1.),
                   ctps_center[3][2]},
        PointADT3D{0.5 * (1. - thickness_up),
                   one_third * (ctps_center[2][1] + 2.),
                   0.5 * (1. - thickness_up)},
        PointADT3D{0.5 * (1. + thickness_up),
                   one_third * (ctps_center[2][1] + 2.),
                   0.5 * (1. - thickness_up)},
        PointADT3D{0.5 * (1. - thickness_up), ADT(1., number_of_derivatives),
                   0.5 * (1. - thickness_up)},
        PointADT3D{0.5 * (1. + thickness_up), ADT(1., number_of_derivatives),
                   0.5 * (1. - thickness_up)},
        ctps_center[6],
        ctps_center[7],
        PointADT3D{ctps_center[6][0], one_third * (2. * ctps_center[7][1] + 1.),
                   ctps_center[6][2]},
        PointADT3D{ctps_center[7][0], one_third * (2. * ctps_center[6][1] + 1.),
                   ctps_center[7][2]},
        PointADT3D{0.5 * (1. - thickness_up),
                   one_third * (ctps_center[2][1] + 2.),
                   0.5 * (1. + thickness_up)},
        PointADT3D{0.5 * (1. + thickness_up),
                   one_third * (ctps_center[2][1] + 2.),
                   0.5 * (1. + thickness_up)},
        PointADT3D{0.5 * (1. - thickness_up), ADT(1., number_of_derivatives),
                   0.5 * (1. + thickness_up)},
        PointADT3D{0.5 * (1. + thickness_up), ADT(1., number_of_derivatives),
                   0.5 * (1. + thickness_up)},
    };

    // Left CTPS
    std::vector<PointADT3D> ctps_left{
        PointADT3D{ADT(0., number_of_derivatives), 0.5 * (1. - thickness_left),
                   0.5 * (1. - thickness_left)},
        PointADT3D{one_third * ctps_center[0][0], 0.5 * (1. - thickness_left),
                   0.5 * (1. - thickness_left)},
        PointADT3D{two_thirds * ctps_center[0][0], ctps_center[0][1],
                   ctps_center[0][2]},
        ctps_center[0],
        PointADT3D{ADT(0., number_of_derivatives), 0.5 * (1. + thickness_left),
                   0.5 * (1. - thickness_left)},
        PointADT3D{one_third * ctps_center[2][0], 0.5 * (1. + thickness_left),
                   0.5 * (1. - thickness_left)},
        PointADT3D{two_thirds * ctps_center[2][0], ctps_center[2][1],
                   ctps_center[2][2]},
        ctps_center[2],
        PointADT3D{ADT(0., number_of_derivatives), 0.5 * (1. - thickness_left),
                   0.5 * (1. + thickness_left)},
        PointADT3D{one_third * ctps_center[4][0], 0.5 * (1. - thickness_left),
                   0.5 * (1. + thickness_left)},
        PointADT3D{two_thirds * ctps_center[4][0], ctps_center[4][1],
                   ctps_center[4][2]},
        ctps_center[4],
        PointADT3D{ADT(0., number_of_derivatives), 0.5 * (1. + thickness_left),
                   0.5 * (1. + thickness_left)},
        PointADT3D{one_third * ctps_center[6][0], 0.5 * (1. + thickness_left),
                   0.5 * (1. + thickness_left)},
        PointADT3D{two_thirds * ctps_center[6][0], ctps_center[6][1],
                   ctps_center[6][2]},
        ctps_center[6]};

    // Right CTPS
    std::vector<PointADT3D> ctps_right{
        ctps_center[1],
        PointADT3D{one_third * (1. + 2. * ctps_center[1][0]), ctps_center[1][1],
                   ctps_center[1][2]},
        PointADT3D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. - thickness_right), 0.5 * (1. - thickness_right)},
        PointADT3D{ADT(1., number_of_derivatives), 0.5 * (1. - thickness_right),
                   0.5 * (1. - thickness_right)},
        ctps_center[3],
        PointADT3D{one_third * (1. + 2. * ctps_center[1][0]), ctps_center[3][1],
                   ctps_center[3][2]},
        PointADT3D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. + thickness_right), 0.5 * (1. - thickness_right)},
        PointADT3D{ADT(1., number_of_derivatives), 0.5 * (1. + thickness_right),
                   0.5 * (1. - thickness_right)},
        ctps_center[5],
        PointADT3D{one_third * (1. + 2. * ctps_center[5][0]), ctps_center[5][1],
                   ctps_center[5][2]},
        PointADT3D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. - thickness_right), 0.5 * (1. + thickness_right)},
        PointADT3D{ADT(1., number_of_derivatives), 0.5 * (1. - thickness_right),
                   0.5 * (1. + thickness_right)},
        ctps_center[7],
        PointADT3D{one_third * (1. + 2. * ctps_center[7][0]), ctps_center[7][1],
                   ctps_center[7][2]},
        PointADT3D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. + thickness_right), 0.5 * (1. + thickness_right)},
        PointADT3D{ADT(1., number_of_derivatives), 0.5 * (1. + thickness_right),
                   0.5 * (1. + thickness_right)}};

    // Front CTPS
    std::vector<PointADT3D> ctps_front{
        PointADT3D{0.5 * (1. - thickness_front), 0.5 * (1. - thickness_front),
                   ADT(0., number_of_derivatives)},
        PointADT3D{0.5 * (1. + thickness_front), 0.5 * (1. - thickness_front),
                   ADT(0., number_of_derivatives)},
        PointADT3D{0.5 * (1. - thickness_front), 0.5 * (1. + thickness_front),
                   ADT(0., number_of_derivatives)},
        PointADT3D{0.5 * (1. + thickness_front), 0.5 * (1. + thickness_front),
                   ADT(0., number_of_derivatives)},
        PointADT3D{0.5 * (1. - thickness_front), 0.5 * (1. - thickness_front),
                   one_third * (ctps_center[0][2])},
        PointADT3D{0.5 * (1. + thickness_front), 0.5 * (1. - thickness_front),
                   one_third * (ctps_center[1][2])},
        PointADT3D{0.5 * (1. - thickness_front), 0.5 * (1. + thickness_front),
                   one_third * (ctps_center[2][2])},
        PointADT3D{0.5 * (1. + thickness_front), 0.5 * (1. + thickness_front),
                   one_third * (ctps_center[3][2])},
        PointADT3D{ctps_center[0][0], ctps_center[0][1],
                   two_thirds * (ctps_center[0][2])},
        PointADT3D{ctps_center[1][0], ctps_center[1][1],
                   two_thirds * (ctps_center[1][2])},
        PointADT3D{ctps_center[2][0], ctps_center[2][1],
                   two_thirds * (ctps_center[2][2])},
        PointADT3D{ctps_center[3][0], ctps_center[3][1],
                   two_thirds * (ctps_center[3][2])},
        ctps_center[0],
        ctps_center[1],
        ctps_center[2],
        ctps_center[3]};

    // Back CTPS
    std::vector<PointADT3D> ctps_back{
        ctps_center[4],
        ctps_center[5],
        ctps_center[6],
        ctps_center[7],
        PointADT3D{ctps_center[4][0], ctps_center[4][1],
                   one_third * (1. + 2. * ctps_center[4][2])},
        PointADT3D{ctps_center[5][0], ctps_center[5][1],
                   one_third * (1. + 2. * ctps_center[5][2])},
        PointADT3D{ctps_center[6][0], ctps_center[6][1],
                   one_third * (1. + 2. * ctps_center[6][2])},
        PointADT3D{ctps_center[7][0], ctps_center[7][1],
                   one_third * (1. + 2. * ctps_center[7][2])},
        PointADT3D{0.5 * (1. - thickness_back), 0.5 * (1. - thickness_back),
                   one_third * (2. + 1. * ctps_center[4][2])},
        PointADT3D{0.5 * (1. + thickness_back), 0.5 * (1. - thickness_back),
                   one_third * (2. + 1. * ctps_center[5][2])},
        PointADT3D{0.5 * (1. - thickness_back), 0.5 * (1. + thickness_back),
                   one_third * (2. + 1. * ctps_center[6][2])},
        PointADT3D{0.5 * (1. + thickness_back), 0.5 * (1. + thickness_back),
                   one_third * (2. + 1. * ctps_center[7][2])},
        PointADT3D{0.5 * (1. - thickness_back), 0.5 * (1. - thickness_back),
                   ADT(1., number_of_derivatives)},
        PointADT3D{0.5 * (1. + thickness_back), 0.5 * (1. - thickness_back),
                   ADT(1., number_of_derivatives)},
        PointADT3D{0.5 * (1. - thickness_back), 0.5 * (1. + thickness_back),
                   ADT(1., number_of_derivatives)},
        PointADT3D{0.5 * (1. + thickness_back), 0.5 * (1. + thickness_back),
                   ADT(1., number_of_derivatives)}};

    // Initialize return value (with one additional spline for the value)
    std::vector<PolyBezierGroup> return_value_group(number_of_derivatives + 1);

    // Construct the Microtile as first component of the group
    return_value_group[0] = PolyBezierGroup{
        Bezier{center_degrees, ExtractValuesFromPointCloud(ctps_center)},
        Bezier{vertical_degrees, ExtractValuesFromPointCloud(ctps_low)},
        Bezier{vertical_degrees, ExtractValuesFromPointCloud(ctps_up)},
        Bezier{horizontal_degrees, ExtractValuesFromPointCloud(ctps_left)},
        Bezier{horizontal_degrees, ExtractValuesFromPointCloud(ctps_right)},
        Bezier{pass_z_degrees, ExtractValuesFromPointCloud(ctps_front)},
        Bezier{pass_z_degrees, ExtractValuesFromPointCloud(ctps_back)}};

    for (std::size_t i_deriv{}; i_deriv < number_of_derivatives; i_deriv++) {
      return_value_group[i_deriv + 1ul] = PolyBezierGroup{
          Bezier{center_degrees,
                 ExtractDerivativeFromPointCloud(ctps_center, i_deriv)},
          Bezier{vertical_degrees,
                 ExtractDerivativeFromPointCloud(ctps_low, i_deriv)},
          Bezier{vertical_degrees,
                 ExtractDerivativeFromPointCloud(ctps_up, i_deriv)},
          Bezier{horizontal_degrees,
                 ExtractDerivativeFromPointCloud(ctps_left, i_deriv)},
          Bezier{horizontal_degrees,
                 ExtractDerivativeFromPointCloud(ctps_right, i_deriv)},
          Bezier{pass_z_degrees,
                 ExtractDerivativeFromPointCloud(ctps_front, i_deriv)},
          Bezier{pass_z_degrees,
                 ExtractDerivativeFromPointCloud(ctps_back, i_deriv)}};
    }

    return return_value_group;
  }

  /// Face to be closed
  enum ClosingFace { X_MIN, X_MAX, Y_MIN, Y_MAX, Z_MIN, Z_MAX };

  /// Provides a closing tile (tile that closes space water tight in a certain
  /// direction)
  static constexpr bool HAS_CLOSING_TILE_DEFINITION = false;
  /// Number of Evaluation Points for Clusing Tile
  static constexpr std::size_t kNumberOfEvaluationPointsClosingTile = 6;
  /// Evaluation Points
  static constexpr std::array<Point3D, kNumberOfEvaluationPointsClosingTile>
      kClosingTileEvaluationPoints{
          Point3D{.5, 0., 0.5},  Point3D{1., 0.5, 0.5}, Point3D{.5, 1., 0.5},
          Point3D{0., 0.5, 0.5}, Point3D{0.5, 0.5, 0.}, Point3D{0.5, 0.5, 1.}};

  /// Generator Function for Derivatives
  static std::vector<PolyBezierGroup> GenerateMicrostructureClosingTile(
      const ClosingFace& closing_face,
      const std::array<ADT, kNumberOfEvaluationPointsClosingTile>&
          evaluations) {
    // Check if all evaluations have the same number of derivatives
    const auto number_of_derivatives = evaluations[0].GetNumberOfDerivatives();
    assert(std::all_of(evaluations.begin(), evaluations.end(),
                       [&number_of_derivatives](const ADT& ev_point) -> bool {
                         return (ev_point.GetNumberOfDerivatives() ==
                                 number_of_derivatives);
                       }));

    // Name abreviations (only readability no purpose)
    const auto& thickness_low = evaluations[0];
    const auto& thickness_right = evaluations[1];
    const auto& thickness_up = evaluations[2];
    const auto& thickness_left = evaluations[3];
    const auto& thickness_front = evaluations[2];
    const auto& thickness_back = evaluations[3];

    // Define degrees
    std::array<std::size_t, 3> base_degrees{1, 1, 1};
    std::array<std::size_t, 3> connector_degrees{1, 1, 2};

    // Init ctps vectors
    std::vector<PointADT3D> ctps_connector;
    std::vector<PointADT3D> ctps_base;

    // Provide abriviations for constant numbers
    const auto zero_value = ADT(0., number_of_derivatives);
    const auto one_value = ADT(1., number_of_derivatives);
    const auto one_quarter = ADT(.25, number_of_derivatives);
    const auto half_value = ADT(0.5, number_of_derivatives);
    const auto three_quarters = ADT(0.75, number_of_derivatives);

    // Provide Control Point Vectors for every possible case except for the
    // edges, assuming that only one edge is closed at a time
    switch (closing_face) {
      case ClosingFace::X_MIN:
        ctps_base = std::vector<PointADT3D>{
            PointADT3D{zero_value, zero_value, zero_value},
            PointADT3D{zero_value, zero_value, one_value},
            PointADT3D{zero_value, one_value, zero_value},
            PointADT3D{zero_value, one_value, one_value},
            PointADT3D{half_value, zero_value, zero_value},
            PointADT3D{half_value, zero_value, one_value},
            PointADT3D{half_value, one_value, zero_value},
            PointADT3D{half_value, one_value, one_value}};

        ctps_connector = std::vector<PointADT3D>{
            PointADT3D{half_value, zero_value, zero_value},
            PointADT3D{half_value, zero_value, one_value},
            PointADT3D{half_value, one_value, zero_value},
            PointADT3D{half_value, one_value, one_value},
            PointADT3D{three_quarters, 0.5 * (1. - thickness_right),
                       0.5 * (1. - thickness_right)},
            PointADT3D{three_quarters, 0.5 * (1. - thickness_right),
                       0.5 * (1. + thickness_right)},
            PointADT3D{three_quarters, 0.5 * (1. + thickness_right),
                       0.5 * (1. - thickness_right)},
            PointADT3D{three_quarters, 0.5 * (1. + thickness_right),
                       0.5 * (1. + thickness_right)},
            PointADT3D{one_value, 0.5 * (1. - thickness_right),
                       0.5 * (1. - thickness_right)},
            PointADT3D{one_value, 0.5 * (1. - thickness_right),
                       0.5 * (1. + thickness_right)},
            PointADT3D{one_value, 0.5 * (1. + thickness_right),
                       0.5 * (1. - thickness_right)},
            PointADT3D{one_value, 0.5 * (1. + thickness_right),
                       0.5 * (1. + thickness_right)}};
        break;
      case ClosingFace::X_MAX:
        ctps_base = std::vector<PointADT3D>{
            PointADT3D{one_value, zero_value, zero_value},
            PointADT3D{one_value, one_value, zero_value},
            PointADT3D{one_value, zero_value, one_value},
            PointADT3D{one_value, one_value, one_value},
            PointADT3D{half_value, zero_value, zero_value},
            PointADT3D{half_value, one_value, zero_value},
            PointADT3D{half_value, zero_value, one_value},
            PointADT3D{half_value, one_value, one_value}};

        ctps_connector = std::vector<PointADT3D>{
            PointADT3D{half_value, zero_value, zero_value},
            PointADT3D{half_value, one_value, zero_value},
            PointADT3D{half_value, zero_value, one_value},
            PointADT3D{half_value, one_value, one_value},
            PointADT3D{one_quarter, 0.5 * (1. - thickness_left),
                       0.5 * (1. - thickness_left)},
            PointADT3D{one_quarter, 0.5 * (1. + thickness_left),
                       0.5 * (1. - thickness_left)},
            PointADT3D{one_quarter, 0.5 * (1. - thickness_left),
                       0.5 * (1. + thickness_left)},
            PointADT3D{one_quarter, 0.5 * (1. + thickness_left),
                       0.5 * (1. + thickness_left)},
            PointADT3D{zero_value, 0.5 * (1. - thickness_left),
                       0.5 * (1. - thickness_left)},
            PointADT3D{zero_value, 0.5 * (1. + thickness_left),
                       0.5 * (1. - thickness_left)},
            PointADT3D{zero_value, 0.5 * (1. - thickness_left),
                       0.5 * (1. + thickness_left)},
            PointADT3D{zero_value, 0.5 * (1. + thickness_left),
                       0.5 * (1. + thickness_left)}};
        break;
      default:
        assert(("Unknown Case", false));
    }
    // Initialize return value (with one additional spline for the value)
    std::vector<PolyBezierGroup> return_value_group(number_of_derivatives + 1);

    // Construct the Microtile as first component of the group
    return_value_group[0] = PolyBezierGroup{
        Bezier{connector_degrees, ExtractValuesFromPointCloud(ctps_connector)},
        Bezier{base_degrees, ExtractValuesFromPointCloud(ctps_base)}};

    for (std::size_t i_deriv{}; i_deriv < number_of_derivatives; i_deriv++) {
      return_value_group[i_deriv + 1ul] = PolyBezierGroup{
          Bezier{connector_degrees,
                 ExtractDerivativeFromPointCloud(ctps_connector, i_deriv)},
          Bezier{base_degrees,
                 ExtractDerivativeFromPointCloud(ctps_base, i_deriv)}};
    }

    return return_value_group;
  }
};

#endif  // EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROTILE_HPP
