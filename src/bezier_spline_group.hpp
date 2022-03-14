#ifndef SRC_BEZIER_SPLINE_GROUP_HPP
#define SRC_BEZIER_SPLINE_GROUP_HPP

#include <cassert>
#include <vector>

#include "bezierManipulation/src/bezier_spline.hpp"
#include "bezierManipulation/src/point.hpp"

namespace beziermanipulation {

/*
 * Group of Bezier splines with same parametric dimension and PointType
 *
 * Predefining the PointType and the parametric dimension is somewhat
 * restricting, but avoids overcomplicated CRTP implementation. Using tuples
 * woud restrict the use of this function to compile time descriptions and
 * prevents use of export import routines.
 *
 * Class inherits from std::vector. This facilitates the use of its methods in
 * later applications. No operator[] overload required etc.
 *
 * In the microstructure use case, BezierSplineGroup can be used both ways,
 * either as the Microtile, or as the deformation function patches, that form
 * the
 *
 * @tparam parametric_dimension parametric dimension of spline
 * @tparam PhysicalPointType Physical Mapping space
 * @tparam ScalarType Scalar used for interpolation function
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::Scalar>
class BezierSplineGroup
    : public std::vector<
          BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>> {
 private:
  using IndexingType = std::size_t;
  using SplineBaseType =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;
  using BaseVector = std::vector<
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>;

public:
  /// Default constructor (to profit from std::vectors implementations)
  constexpr BezierSplineGroup() = default;

  /// Initializer list overload
  template <typename... Splines>
  constexpr BezierSplineGroup(const Splines &... splines)
      : BaseVector{splines...} {}


  /// Check if group fits unit cube
  constexpr bool fits_unit_cube() const;

  /// Maximum
  constexpr PhysicalPointType maximum() const;

  /// Minimum corner
  constexpr PhysicalPointType minimum() const;

  /// Fit to unit_cube
  constexpr BezierSplineGroup &fit_to_unit_cube();

  /// Compose with single Spline
  constexpr BezierSplineGroup compose(
      const SplineBaseType &inner_function) const;
  /// Compose with Splinegroup
  constexpr BezierSplineGroup compose(
      const BezierSplineGroup &inner_function_group) const;
  // @todo Overload operators if required
};

#include "bezierManipulation/src/bezier_spline_group.inc"

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_SPLINE_GROUP_HPP
