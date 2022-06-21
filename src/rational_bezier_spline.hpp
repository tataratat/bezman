#ifndef SRC_RATIONAL_BEZIER_SPLINE_HPP
#define SRC_RATIONAL_BEZIER_SPLINE_HPP

#include "bezman/src/bezier_spline.hpp"

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class RationalBezierSplineGroup;

/*
 * Class describing rational Bezier-Splines
 *
 * @tparam parametric_dimenasion parametric dimension of the spline
 * @tparam PhysicalPointType Type of the control points that are used for
 * interpolations
 * @tparam ScalarType default scalar type used to in the physical domain
 * (e.g. double / AD-Type)
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::Scalar>
class RationalBezierSpline {
 private:
  // Aliases
  using IndexingType = std::size_t;
  using PointTypePhysical_ = PhysicalPointType;
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;

  /*!
   * Numerator Spline, with control points multiplied with their respective
   * weights
   */
  BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
      weighted_spline_;

  /*!
   * Denominator Spline contains the weight function
   */
  BezierSpline<parametric_dimension, ScalarType, ScalarType> weight_function_;

  /**
   * @brief  Wrapper function that checks if the Numerator and Denominator
   * splines are still compatible after operations are performed
   */
  constexpr bool CheckSplineCompatibility() const {
    bool compatible =
        (weighted_spline_.GetDegrees() == weight_function_.GetDegrees());
    return compatible;
  }

 public:
  /// Make Parametric dimension publicly available
  static constexpr IndexingType kParametricDimensions = parametric_dimension;

  /// Copy constructor
  constexpr RationalBezierSpline(const RationalBezierSpline& bezier_spline) =
      default;

  /// Transform Bezier into rational Bezier
  template <
      typename BezierSplineT,
      std::enable_if_t<type_traits::isBezier_v<OutputType>, void*> = nullptr>
  constexpr RationalBezierSpline(const BezierSpline& bezier_spline)
      : weighted_spline_{bezier_spline} {
    weight_function_(
        // Degrees can be taken from weighted spline
        weighted_spline_.GetDegrees(),
        // Generate control weights (just 1 for every point)
        std::move(std::vector(
            // Size of control-point vector
            weighted_spline.NumberOfControlPoints,
            // Value
            static_cast<ScalarType>(1.))));
    assert(CheckSplineCompatibility());
  }

  /// Empty constructor
  constexpr RationalBezierSpline() = default;

  /// Move operator
  constexpr RationalBezierSpline& operator=(RationalBezierSpline&& rhs) =
      default;

  /// Move operator
  constexpr RationalBezierSpline& operator=(const RationalBezierSpline& rhs) =
      default;

  /// Getter for Degrees
  constexpr const std::array<std::size_t, parametric_dimension>& GetDegrees()
      const {
    return weighted_spline.GetDegrees;
  }

  // @NOTE These functions all pass arguments and perform very limited
  // operations

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& ControlPoint(const T... index);

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& WeightedControlPoint(
      const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& WeightedControlPoint(const T... index);

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& Weight(const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& Weight(const T... index);

  // @NOTE Theoretically, these can all be implemented for arrays as well

  //-------------------

  // @NOTE Essential operations

  /// Order elevation along a specific parametric dimension
  constexpr RationalBezierSpline& OrderElevateAlongParametricDimension(
      const IndexingType par_dim);

  // @Note as far as I can see, there is no simple way to perform this
  // operation and we will have to resort to quotient-rule derivatives
  /// Derivative along a specific parametric dimension
  constexpr RationalBezierSpline DerivativeWRTParametricDimension(
      const IndexingType par_dim) const;

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ Evaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using the de Casteljau algorithm
  template <typename... T>
  constexpr PointTypePhysical_ Evaluate(const T&... par_coords) const {
    return Evaluate(PointTypeParametric_{par_coords...});
  }
  /// Evaluate the spline via the explicit precomputation of bernstein values
  constexpr PointTypePhysical_ ForwardEvaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using explicit precomputation of bernstein values
  template <typename... T>
  constexpr PointTypePhysical_ ForwardEvaluate(const T&... par_coords) const {
    return ForwardEvaluate(PointTypeParametric_{par_coords...});
  }

  /*
   * Functional Composition between two splines
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and the
   * function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension_inner_spline,
                                 PhysicalPointType,
                                 decltype(ScalarType_{} * ScalarRHS{})>
  Compose(const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                             ScalarRHS>& inner_function) const;

  /*
   * Split the Bezier Spline into two distinct subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over two splines.
   */
  constexpr RationalBezierSplineGroup<parametric_dimension, PhysicalPointType,
                                      ScalarType>
  SplitAtPosition(const ScalarType& splitting_plane,
                  const IndexingType splitting_dimension = 0) const;
}
#endif  // SRC_RATIONAL_BEZIER_SPLINE_HPP
