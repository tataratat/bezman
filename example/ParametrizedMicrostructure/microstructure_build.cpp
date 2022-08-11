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

/*
 * This Example illustrates, how the BezierManipulation suite can be used to
 * build parametrized microstructures. It is also used to show the capabilities
 * of deriving the structure with respect to the design variables, using
 * algorithmic differentiation.
 *
 * Implementationwise, the code is structured as follows. There are three
 * classes describing the geometry and the parametrization. The MicortileExample
 * class, describes the parametrized Inner-spline, that is later inserted into
 * the structure. It contains a function, that creates a microtile, based on a
 * set of parameters. As the parameters are based on a function evaluation, it
 * also provides the points, where this funciton is to be evaluated (in the case
 * of a crosstile, these are the tips of the crosstile-branches). The
 * DeformationFunctionExample class provides the outer geometry. As it is formed
 * through a group of Bezier splines, it also holds segmentation information,
 * i.e., how many splines are in each parametric dimension. This could also be
 * generalized with a function, returning local coordinates based on the current
 * index. The ValueFunction class finally provides the parametrization function,
 * that returns the local Parameters, based on the position of the spline within
 * the structure and a set of superordinate parameters. These three functions
 * are then inserted as Template Arguments in the Microstructure-Generator
 * class, which performes the compositions and derivations.
 *
 * The results are exported in xml format and can be plotted using the python
 * scripts with gustav
 *
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include "bezman/src/bezier_group.hpp"
#include "bezman/src/utils/computational_differentiation/algo_diff_type.hpp"
#include "bezman/src/utils/export.hpp"

using namespace bezman;

// Define aliases to abriviate code
using ADT = utils::computational_differentiation::AlgoDiffType<double>;
using PointADT2D = Point<2, ADT>;
using Point2D = Point<2, double>;
using Bezier = BezierSpline<2, Point2D, double>;
using PolynomialBezierGroup = BezierGroup<Bezier>;

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
class MicrotileExample {
 private:
  /**
   * @brief Extracts the Derivatives from the derivative vector of each control
   * point and defines new sets of splines, which  are then used to determine
   * the derivatives of the microtiles
   */
  static std::vector<Point2D> ExtractDerivativeFromPointCloud(
      const std::vector<PointADT2D>& points, const std::size_t& derivative) {
    // Initialize return type
    std::vector<Point2D> return_vector(points.size());

    for (std::size_t i_point{}; i_point < points.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < PointADT2D::kSpatialDimension;
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
   * @return std::vector<Point2D>
   */
  static std::vector<Point2D> ExtractValuesFromPointCloud(
      const std::vector<PointADT2D>& points) {
    // Initialize return type
    std::vector<Point2D> return_vector(points.size());
    for (std::size_t i_point{}; i_point < points.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < 2; i_dim++) {
        return_vector[i_point][i_dim] = points[i_point][i_dim].GetValue();
      }
    }
    return return_vector;
  }

 public:
  /// Number of Evaluation Points
  static constexpr std::size_t kNumberOfEvaluationPoints = 4;
  /// Evaluation Points
  static constexpr std::array<Point2D, kNumberOfEvaluationPoints>
      kEvaluationPoints{Point2D{.5, 0.}, Point2D{1., 0.5}, Point2D{.5, 1.},
                        Point2D{0., 0.5}};

  /// Generator Function for Derivatives
  static std::vector<PolynomialBezierGroup> GenerateMicrostructureDerivatives(
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

    // Center dimensions
    const auto center_width = 0.5 * (evaluations[0] + evaluations[2]);
    const auto center_height = 0.5 * (evaluations[3] + evaluations[1]);

    // Name abreviations (only readability no purpose)
    const auto& thickness_low = evaluations[0];
    const auto& thickness_right = evaluations[1];
    const auto& thickness_up = evaluations[2];
    const auto& thickness_left = evaluations[3];

    // Define degrees
    std::array<std::size_t, 2> center_degrees{1, 1};
    std::array<std::size_t, 2> horizontal_degrees{3, 1};
    std::array<std::size_t, 2> vertical_degrees{1, 3};

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
    std::vector<PointADT2D> ctps_center{
        PointADT2D{0.5 * (1. - center_width), 0.5 * (1. - center_height)},
        PointADT2D{0.5 * (1. + center_width), 0.5 * (1. - center_height)},
        PointADT2D{0.5 * (1. - center_width), 0.5 * (1. + center_height)},
        PointADT2D{0.5 * (1. + center_width), 0.5 * (1. + center_height)}};

    // Lower CTPS
    std::vector<PointADT2D> ctps_low{
        PointADT2D{0.5 * (1. - thickness_low), ADT(0., number_of_derivatives)},
        PointADT2D{0.5 * (1. + thickness_low), ADT(0., number_of_derivatives)},
        PointADT2D{0.5 * (1. - thickness_low), one_third * ctps_center[0][1]},
        PointADT2D{0.5 * (1. + thickness_low), one_third * ctps_center[0][1]},
        PointADT2D{ctps_center[0][0], two_thirds * ctps_center[0][1]},
        PointADT2D{ctps_center[1][0], two_thirds * ctps_center[0][1]},
        ctps_center[0],
        ctps_center[1]};

    // Upper CTPS
    std::vector<PointADT2D> ctps_up{
        ctps_center[2],
        ctps_center[3],
        PointADT2D{ctps_center[2][0],
                   one_third * (2. * ctps_center[2][1] + 1.)},
        PointADT2D{ctps_center[3][0],
                   one_third * (2. * ctps_center[2][1] + 1.)},
        PointADT2D{0.5 * (1. - thickness_up),
                   one_third * (ctps_center[2][1] + 2.)},
        PointADT2D{0.5 * (1. + thickness_up),
                   one_third * (ctps_center[2][1] + 2.)},
        PointADT2D{0.5 * (1. - thickness_up), ADT(1., number_of_derivatives)},
        PointADT2D{0.5 * (1. + thickness_up), ADT(1., number_of_derivatives)}};

    // Left CTPS
    std::vector<PointADT2D> ctps_left{
        PointADT2D{ADT(0., number_of_derivatives), 0.5 * (1. - thickness_left)},
        PointADT2D{one_third * ctps_center[0][0], 0.5 * (1. - thickness_left)},
        PointADT2D{two_thirds * ctps_center[0][0], ctps_center[0][1]},
        ctps_center[0],
        PointADT2D{ADT(0., number_of_derivatives), 0.5 * (1. + thickness_left)},
        PointADT2D{one_third * ctps_center[2][0], 0.5 * (1. + thickness_left)},
        PointADT2D{two_thirds * ctps_center[2][0], ctps_center[2][1]},
        ctps_center[2]};

    // Right CTPS
    std::vector<PointADT2D> ctps_right{
        ctps_center[1],
        PointADT2D{one_third * (1. + 2. * ctps_center[1][0]),
                   ctps_center[1][1]},
        PointADT2D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. - thickness_right)},
        PointADT2D{ADT(1., number_of_derivatives),
                   0.5 * (1. - thickness_right)},
        ctps_center[3],
        PointADT2D{one_third * (1. + 2. * ctps_center[1][0]),
                   ctps_center[3][1]},
        PointADT2D{one_third * (2. + ctps_center[1][0]),
                   0.5 * (1. + thickness_right)},
        PointADT2D{ADT(1., number_of_derivatives),
                   0.5 * (1. + thickness_right)}};

    // Initialize return value (with one additional spline for the value)
    std::vector<PolynomialBezierGroup> return_value_group(
        number_of_derivatives + 1);

    // Construct the Microtile as first component of the group
    return_value_group[0] = PolynomialBezierGroup{
        Bezier{center_degrees, ExtractValuesFromPointCloud(ctps_center)},
        Bezier{vertical_degrees, ExtractValuesFromPointCloud(ctps_low)},
        Bezier{vertical_degrees, ExtractValuesFromPointCloud(ctps_up)},
        Bezier{horizontal_degrees, ExtractValuesFromPointCloud(ctps_left)},
        Bezier{horizontal_degrees, ExtractValuesFromPointCloud(ctps_right)}};

    for (std::size_t i_deriv{}; i_deriv < number_of_derivatives; i_deriv++) {
      return_value_group[i_deriv + 1ul] = PolynomialBezierGroup{
          Bezier{center_degrees,
                 ExtractDerivativeFromPointCloud(ctps_center, i_deriv)},
          Bezier{vertical_degrees,
                 ExtractDerivativeFromPointCloud(ctps_low, i_deriv)},
          Bezier{vertical_degrees,
                 ExtractDerivativeFromPointCloud(ctps_up, i_deriv)},
          Bezier{horizontal_degrees,
                 ExtractDerivativeFromPointCloud(ctps_left, i_deriv)},
          Bezier{horizontal_degrees,
                 ExtractDerivativeFromPointCloud(ctps_right, i_deriv)}};
    }

    return return_value_group;
  }

  /// Face to be closed
  enum ClosingFace { X_MIN, X_MAX, Y_MIN, Y_MAX };

  /// Provides a closing tile (tile that closes space water tight in a certain
  /// direction)
  static constexpr bool HAS_CLOSING_TILE_DEFINITION = true;
  /// Number of Evaluation Points for Clusing Tile
  static constexpr std::size_t kNumberOfEvaluationPointsClosingTile = 4;
  /// Evaluation Points
  static constexpr std::array<Point2D, kNumberOfEvaluationPointsClosingTile>
      kClosingTileEvaluationPoints{Point2D{.5, 0.}, Point2D{1., 0.5},
                                   Point2D{.5, 1.}, Point2D{0., 0.5}};

  /// Generator Function for Derivatives
  static std::vector<PolynomialBezierGroup> GenerateMicrostructureClosingTile(
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

    // Define degrees
    std::array<std::size_t, 2> base_degrees{1, 1};
    std::array<std::size_t, 2> connector_degrees{1, 2};

    std::vector<PointADT2D> ctps_connector;
    std::vector<PointADT2D> ctps_base, ctps_base_lhs, ctps_base_rhs;

    // Provide Control Point Vectors for every possible case except for the
    // edges, assuming that only one edge is closed at a time
    switch (closing_face) {
      case ClosingFace::Y_MIN:
        ctps_connector = std::vector<PointADT2D>{
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{0.5 * (1. - thickness_up),
                       ADT(0.75, number_of_derivatives)},
            PointADT2D{0.5 * (1. + thickness_up),
                       ADT(0.75, number_of_derivatives)},
            PointADT2D{0.5 * (1. - thickness_up),
                       ADT(1., number_of_derivatives)},
            PointADT2D{0.5 * (1. + thickness_up),
                       ADT(1., number_of_derivatives)}};

        ctps_base = std::vector<PointADT2D>{
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)}};

        ctps_base_lhs = std::vector<PointADT2D>{
            PointADT2D{ADT(0.0, number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(0.0, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)}};

        ctps_base_rhs = std::vector<PointADT2D>{
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(1., number_of_derivatives),
                       ADT(.0, number_of_derivatives)},
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{ADT(1., number_of_derivatives),
                       ADT(0.5, number_of_derivatives)}};
        break;
      case ClosingFace::Y_MAX:
        ctps_connector = std::vector<PointADT2D>{
            PointADT2D{0.5 * (1. - thickness_low),
                       ADT(0., number_of_derivatives)},
            PointADT2D{0.5 * (1. + thickness_low),
                       ADT(0., number_of_derivatives)},
            PointADT2D{0.5 * (1. - thickness_low),
                       ADT(0.25, number_of_derivatives)},
            PointADT2D{0.5 * (1. + thickness_low),
                       ADT(0.25, number_of_derivatives)},
            PointADT2D{ADT(0.05, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)},
            PointADT2D{ADT(0.95, number_of_derivatives),
                       ADT(0.5, number_of_derivatives)}};

        ctps_base =
            std::vector<PointADT2D>{PointADT2D{ADT(0.05, number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(0.95, number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(0.05, number_of_derivatives),
                                               ADT(1., number_of_derivatives)},
                                    PointADT2D{ADT(0.95, number_of_derivatives),
                                               ADT(1., number_of_derivatives)}};

        ctps_base_lhs =
            std::vector<PointADT2D>{PointADT2D{ADT(0., number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(0.05, number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(0., number_of_derivatives),
                                               ADT(1., number_of_derivatives)},
                                    PointADT2D{ADT(0.05, number_of_derivatives),
                                               ADT(1., number_of_derivatives)}};

        ctps_base_rhs =
            std::vector<PointADT2D>{PointADT2D{ADT(0.95, number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(1., number_of_derivatives),
                                               ADT(0.5, number_of_derivatives)},
                                    PointADT2D{ADT(0.95, number_of_derivatives),
                                               ADT(1., number_of_derivatives)},
                                    PointADT2D{ADT(1., number_of_derivatives),
                                               ADT(1., number_of_derivatives)}};
        break;
      default:
        assert(("Unknown Case", false));
    }
    // Initialize return value (with one additional spline for the value)
    std::vector<PolynomialBezierGroup> return_value_group(
        number_of_derivatives + 1);

    // Construct the Microtile as first component of the group
    return_value_group[0] = PolynomialBezierGroup{
        Bezier{connector_degrees, ExtractValuesFromPointCloud(ctps_connector)},
        Bezier{base_degrees, ExtractValuesFromPointCloud(ctps_base)},
        Bezier{base_degrees, ExtractValuesFromPointCloud(ctps_base_lhs)},
        Bezier{base_degrees, ExtractValuesFromPointCloud(ctps_base_rhs)}};

    for (std::size_t i_deriv{}; i_deriv < number_of_derivatives; i_deriv++) {
      return_value_group[i_deriv + 1ul] = PolynomialBezierGroup{
          Bezier{connector_degrees,
                 ExtractDerivativeFromPointCloud(ctps_connector, i_deriv)},
          Bezier{base_degrees,
                 ExtractDerivativeFromPointCloud(ctps_base, i_deriv)},
          Bezier{base_degrees,
                 ExtractDerivativeFromPointCloud(ctps_base_lhs, i_deriv)},
          Bezier{base_degrees,
                 ExtractDerivativeFromPointCloud(ctps_base_rhs, i_deriv)}};
    }

    return return_value_group;
  }
};

/**
 * @brief Deformation function defining the outer contour of a spline
 *
 * Approximates a circle with a given number of segments in each parametric
 * dimension. Circles are approximated with second order B-Spline segments.
 *
 */
class DeformationFunctionExample {
 private:
  /// Precalculated values
  constexpr static const double PI = std::acos(-1.);

 public:
  // Quarter Circle Dimensions
  //  static constexpr const double innerR{8.}, outerR{20.}, arc_degrees{2 *
  //  PI};
  static constexpr const double innerR{8.}, outerR{20.}, arc_degrees{0.5 * PI};

  /// Number of segments in each parametric dimension
  static constexpr const int kNumberOfXSegments{5};
  static constexpr const int kNumberOfYSegments{5};

  /// For external access
  static constexpr const std::array<int, 2> kSegmentsPerParametricDimension{
      kNumberOfXSegments, kNumberOfYSegments};

  /**
   * @brief Create Circle Deformation function with given number of segments
   * per parametric dimension
   */
  static PolynomialBezierGroup Create() {
    // split planes to segment circle in radial direction (even samples)
    std::vector<double> y_knot_lines(kNumberOfYSegments - 1);
    const double y_seg_length = 1. / static_cast<double>(kNumberOfYSegments);
    for (unsigned int i_y_segment{1}; i_y_segment < kNumberOfYSegments;
         i_y_segment++) {
      y_knot_lines[i_y_segment - 1] = y_seg_length * i_y_segment;
    }

    // Initialize return value
    PolynomialBezierGroup ringsegments{kNumberOfXSegments * kNumberOfYSegments};

    // Precompute values that are required multiple times
    const std::array<std::size_t, 2> degrees{2, 1};
    const double degrees_per_segment = arc_degrees / kNumberOfXSegments;
    const double excentricity_of_middle_points =
        1. / std::sin(PI / 2. - degrees_per_segment / 2.);

    for (int i_x_segment{}; i_x_segment < kNumberOfXSegments; i_x_segment++) {
      const double startsin = std::sin(degrees_per_segment * i_x_segment);
      const double startcos = std::cos(degrees_per_segment * i_x_segment);
      const double middlesin =
          std::sin(degrees_per_segment * (i_x_segment + .5));
      const double middlecos =
          std::cos(degrees_per_segment * (i_x_segment + .5));
      const double endsin = std::sin(degrees_per_segment * (i_x_segment + 1.));
      const double endcos = std::cos(degrees_per_segment * (i_x_segment + 1.));

      // Define the control points
      std::vector<Point2D> ctps{
          Point2D{startcos * outerR, startsin * outerR},
          Point2D{middlecos * excentricity_of_middle_points * outerR,
                  middlesin * excentricity_of_middle_points * outerR},
          Point2D{endcos * outerR, endsin * outerR},
          Point2D{startcos * innerR, startsin * innerR},
          Point2D{middlecos * excentricity_of_middle_points * innerR,
                  middlesin * excentricity_of_middle_points * innerR},
          Point2D{endcos * innerR, endsin * innerR}};
      // Define and split up circle segment
      const auto CircleWedge =
          Bezier{degrees, ctps}.SplitAtPosition(y_knot_lines, 1);
      // Assign splines to Splinegroup
      for (int i_y_segment{}; i_y_segment < kNumberOfYSegments; i_y_segment++) {
        ringsegments[i_x_segment + kNumberOfXSegments * i_y_segment] =
            CircleWedge[i_y_segment];
      }
    }
    return ringsegments;
  }
};

/**
 * @brief Value Field Function defined in the parametric domain of the
 * deformation function
 *
 * This function is based on Algorithmic differentiation types, which calculate
 * the derivative of the value field function at the same time.
 *
 */
class ValueFieldExample {
 public:
  // Number of super Parameters
  static constexpr const std::size_t kNumberOfSuperParameters{6};

  // Some super parameters
  const std::array<ADT, kNumberOfSuperParameters> kSuperControlPoints{
      ADT(0.3, kNumberOfSuperParameters, 0),
      ADT(0.4, kNumberOfSuperParameters, 1),
      ADT(0.3, kNumberOfSuperParameters, 2),
      ADT(0.25, kNumberOfSuperParameters, 3),
      ADT(0.25, kNumberOfSuperParameters, 4),
      ADT(0.25, kNumberOfSuperParameters, 5)};

  /**
   * @brief Evaluate the Value field function
   *
   * Here, the value field function are the basis functions of a B-Spline
   * interpolating between the SuperControlPoints. The B-Spline has a knot in
   * the center, meaning that the basis functions have only local support.
   */
  template <std::size_t array_size>
  std::array<ADT, array_size> Evaluate(
      const std::array<Point2D, array_size>& evaluation_points) const {
    // Initialize return value
    std::array<ADT, array_size> evaluations{};
    // Assign based on simple linear B-Spline basis
    for (std::size_t i_point{}; i_point < array_size; i_point++) {
      const double &x{evaluation_points[i_point][0]},
          y{evaluation_points[i_point][1]};
      evaluations[i_point] =
          (x < 0.5 ? (1 - 2. * x) : 0.) * (1 - y) * kSuperControlPoints[0] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * (1 - y) *
              kSuperControlPoints[1] +
          (x < 0.5 ? 0. : 2 * x - 1.) * (1 - y) * kSuperControlPoints[2] +
          (x < 0.5 ? (1 - 2. * x) : 0.) * y * kSuperControlPoints[3] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * y * kSuperControlPoints[4] +
          (x < 0.5 ? 0. : 2 * x - 1) * y * kSuperControlPoints[5];
    }
    return evaluations;
  }
};
/**
 * @brief Generator of the parametrized microstructure
 *
 * Constructs the parametrized microstructure along with its derivatives
 *
 * @tparam Microtile            Microtile Example (Parametrized function)
 * @tparam DeformationFunction  Deformation function as a PolynomialBezierGroup
 * @tparam ValueField           Value function that is used to determine the
 * local Microtileparameters
 */
template <typename Microtile, typename DeformationFunction, typename ValueField>
class ParametrizedComposition {
 private:
  // TransformPoints_ transforms the points provided by
  // Microtile::kEvaluationPoints these pints are transformed by a linear
  // interpolation funciton and set into the Valuefield function to retrieve the
  // local parameters  of the microtile function.
  std::array<Point2D, Microtile::kNumberOfEvaluationPoints> TransformPoints_(
      const Point2D& cornermin, const Point2D& cornermax,
      const std::array<Point2D, Microtile::kNumberOfEvaluationPoints>&
          evaluation_points) const {
    std::array<Point2D, Microtile::kNumberOfEvaluationPoints> return_value{};
    const auto difference = cornermax - cornermin;
    for (std::size_t i{}; i < Microtile::kNumberOfEvaluationPoints; i++) {
      return_value[i] =
          cornermin + Point2D{evaluation_points[i][0] * difference[0],
                              evaluation_points[i][1] * difference[1]};
    }
    return return_value;
  }

  // Instance of the value field function, this is required as the vector class
  // can not be used at compile time in the cxx17 standard. This feature is
  // added in cxx20. (For this reason I am thinking of changing the standard in
  // the long run)
  ValueField value_field{};

 public:
  // Compose microstructure
  std::vector<PolynomialBezierGroup> ComposeMicrostructureAndDerivatives()
      const {
    // Loop over the externnal splines in the group
    const PolynomialBezierGroup deformation_function =
        DeformationFunction::Create();
    utils::Export::GuessByExtension(deformation_function,
                                    "deformation_function.xml");

    // Initialize return value
    std::vector<PolynomialBezierGroup> return_value(
        ValueFieldExample::kNumberOfSuperParameters + 1u);

    // Every spline in the mirostructure defines a specific range of the
    // total spline group. We consider an even splitting, meaning that that
    // the parametric sections have the same size. As a first step, the
    // current field needs to be identified.
    const double dx =
        1. / static_cast<double>(
                 DeformationFunction::kSegmentsPerParametricDimension[0]);
    const double dy =
        1. / static_cast<double>(
                 DeformationFunction::kSegmentsPerParametricDimension[1]);
    for (std::size_t i_def_x{};
         i_def_x < DeformationFunction::kSegmentsPerParametricDimension[0];
         i_def_x++) {
      for (std::size_t i_def_y{};
           i_def_y < DeformationFunction::kSegmentsPerParametricDimension[1];
           i_def_y++) {
        std::vector<PolynomialBezierGroup> microtile_vector;
        if constexpr (Microtile::HAS_CLOSING_TILE_DEFINITION) {
          // If a closing tile is defined we can create watertight structures
          if (i_def_y == 0) {
            const auto transformed_points = TransformPoints_(
                Point2D{i_def_x * dx, i_def_y * dy},
                Point2D{(i_def_x + 1) * dx, (i_def_y + 1) * dy},
                Microtile::kClosingTileEvaluationPoints);
            // Evaluate microtile based on the generator provided by Microtile
            microtile_vector = Microtile::GenerateMicrostructureClosingTile(
                Microtile::ClosingFace::Y_MIN,
                value_field.Evaluate(transformed_points));
          } else if (i_def_y ==
                     DeformationFunction::kSegmentsPerParametricDimension[1] -
                         1) {
            const auto transformed_points = TransformPoints_(
                Point2D{i_def_x * dx, i_def_y * dy},
                Point2D{(i_def_x + 1) * dx, (i_def_y + 1) * dy},
                Microtile::kClosingTileEvaluationPoints);
            // Evaluate microtile based on the generator provided by Microtile
            microtile_vector = Microtile::GenerateMicrostructureClosingTile(
                Microtile::ClosingFace::Y_MAX,
                value_field.Evaluate(transformed_points));
          } else {
            // If no closing tile is defined, then we set normal tiles in the
            // boundary elements
            const auto transformed_points = TransformPoints_(
                Point2D{i_def_x * dx, i_def_y * dy},
                Point2D{(i_def_x + 1) * dx, (i_def_y + 1) * dy},
                Microtile::kEvaluationPoints);
            // Evaluate microtile based on the generator provided by Microtile
            microtile_vector = Microtile::GenerateMicrostructureDerivatives(
                value_field.Evaluate(transformed_points));
          }
        } else {
          // If no closing tile is defined, then we set normal tiles in the
          // boundary elements
          const auto transformed_points =
              TransformPoints_(Point2D{i_def_x * dx, i_def_y * dy},
                               Point2D{(i_def_x + 1) * dx, (i_def_y + 1) * dy},
                               Microtile::kEvaluationPoints);
          // Evaluate microtile based on the generator provided by Microtile
          microtile_vector = Microtile::GenerateMicrostructureDerivatives(
              value_field.Evaluate(transformed_points));
        }

        // Compose the microstructure in first position of vector
        const std::size_t combined_index =
            i_def_x +
            DeformationFunction::kSegmentsPerParametricDimension[0] * i_def_y;

        // Construct the actual Composition (microstructure)
        return_value[0] +=
            deformation_function[combined_index].Compose(microtile_vector[0]);

        /// Construct the corresponding derivatives
        // Precompute the derivatives of the deformation function as they are
        // required multiple times (for every internal control points of the
        // value field funciton)
        // Here it is hard-coded for the two dimensional case
        auto microstructure_derivative_first_par_dim =
            deformation_function[combined_index]
                .DerivativeWRTParametricDimension(0)
                .Compose(microtile_vector[0]);
        auto microstructure_derivative_second_par_dim =
            deformation_function[combined_index]
                .DerivativeWRTParametricDimension(1)
                .Compose(microtile_vector[0]);
        for (std::size_t i_derivative{};
             i_derivative < value_field.kNumberOfSuperParameters;
             i_derivative++) {
          return_value[i_derivative + 1] +=
              (microstructure_derivative_first_par_dim.MultiplyComponentwise(
                   microtile_vector[i_derivative].ExtractDimension(0)))
                  .AddComponentwise(
                      microstructure_derivative_second_par_dim
                          .MultiplyComponentwise(
                              microtile_vector[i_derivative].ExtractDimension(
                                  1)));
        }
      }
    }
    return return_value;
  }

  // Compose Derivatives
};

int main() {
  const auto micro_structure_generator =
      ParametrizedComposition<MicrotileExample, DeformationFunctionExample,
                              ValueFieldExample>{};

  const auto test_composition =
      micro_structure_generator.ComposeMicrostructureAndDerivatives();

  utils::Export::GuessByExtension(test_composition[0],
                                  "composed_microstructure.json");
  utils::Export::GuessByExtension(test_composition[0],
                                  "composed_microstructure.mesh");
  ///   const auto microstructure =
  ///   MicrotileExample::GenerateMicrostructureDerivatives(
  ///       std::array<ADT,4> {
  ///         ADT(0.1,1,0),
  ///         ADT(0.1,1,0),
  ///         ADT(0.1,1,0),
  ///         ADT(0.1,1,0)}
  ///       );
  ///
  ///   utils::Export::GuessByExtension(microstructure[0], "microtile.json");
  ///
  ///   for (std::size_t i_deriv{};
  ///        i_deriv < ValueFieldExample::kNumberOfSuperParameters; i_deriv++) {
  ///     utils::Export::GuessByExtension(test_composition[i_deriv + 1],
  ///                                     std::string("composed_microstructure_")
  ///                                     +
  ///                                         std::to_string(i_deriv) +
  ///                                         std::string(".json"));
  ///   }
  return 0;
}
