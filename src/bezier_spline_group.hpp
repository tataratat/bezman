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
  constexpr bool FitsIntoUnitCube() const;

  /// Maximum
  constexpr PhysicalPointType MaximumCorner() const;

  /// Minimum corner
  constexpr PhysicalPointType MinimumCorner() const;

  /// Fit to unit_cube
  constexpr BezierSplineGroup &FitToUnitCube();

  /// Compose with single Spline
  constexpr BezierSplineGroup Compose(
      const SplineBaseType &inner_function) const;

  /// Compose with Splinegroup
  constexpr BezierSplineGroup Compose(
      const BezierSplineGroup &inner_function_group) const;

  // + Operator for concatenation
  constexpr BezierSplineGroup operator+(const BezierSplineGroup &rhs) const {
    BezierSplineGroup combination{(*this)};
    combination += rhs;
    return combination;
  }

  // + Operator for concatenation
  constexpr BezierSplineGroup &operator+=(const BezierSplineGroup &rhs) {
    (*this).insert(this->begin(), rhs.begin(), rhs.end());
    return (*this);
  }

  // + Operator for concatenation
  constexpr BezierSplineGroup operator+(
      const PhysicalPointType &translation) const {
    BezierSplineGroup combination{(*this)};
    combination += translation;
    return combination;
  }

  // + Operator for concatenation
  constexpr BezierSplineGroup &operator+=(
      const PhysicalPointType &translation) {
    for (IndexingType i_spline{}; i_spline < this->size(); i_spline++) {
      (*this)[i_spline] += translation;
    }
    return (*this);
  }
};

#include "bezierManipulation/src/bezier_spline_group.inc"

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_SPLINE_GROUP_HPP
