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

#ifndef SRC_BEZIER_SPLINE_HPP
#define SRC_BEZIER_SPLINE_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <vector>

#include "bezman/src/point.hpp"
#include "bezman/src/utils/algorithms/bernstein_polynomial.hpp"
#include "bezman/src/utils/algorithms/recursive_combine.hpp"
#include "bezman/src/utils/fastbinomialcoefficient.hpp"
#include "bezman/src/utils/logger.hpp"
#include "bezman/src/utils/type_traits/is_bezier_spline.hpp"
#include "bezman/src/utils/type_traits/is_rational_bezier_spline.hpp"

namespace bezman {

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class RationalBezierSpline;

/*
 * Class describing BezierSplines
 *
 * Core of the library
 *
 * @tparam parametric_dimenasion parametric dimension of the spline
 * @tparam PhysicalPointType Type of the control points that are used for
 * interpolations
 * @tparam ScalarType default scalar type used to in the physical domain
 * (e.g. double / AD-Type)
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::Scalar>
class BezierSpline {
  // Friend declaration
  template <std::size_t parent_parametric_dimension,
            typename ParentPhysicalPointType, typename ParentScalarType>
  friend class RationalBezierSpline;

 private:
  using IndexingType = std::size_t;
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;

  constexpr void UpdateIndexOffsets_();

  /// Polynomial degrees
  std::array<IndexingType, parametric_dimension> degrees{};
  /// Number of control points
  IndexingType number_of_control_points{};

  /**
   * @brief Templated Fallback function for ComposeNumeratorSpline and
   * ComposeNumeratorSensitivity as they share most or their code
   *
   * Check out their respective documentation for more information
   */
  template <bool compute_sensitivities, typename SplineType>
  constexpr auto ComposeNumeratorHelper(const SplineType &inner_function) const;

  /**
   * @brief Compose the Numerator Function of a polynomial and rational spline
   *
   * This function composes the Numerator only, so it can be reused to work with
   * rational-rational-spline-compositions
   */
  template <typename SplineType>
  constexpr auto ComposeNumeratorSpline(
      const SplineType &inner_function) const {
    return ComposeNumeratorHelper<false>(inner_function);
  }

  /**
   * @brief Compose the Sensitivity of the Numerator Function of a polynomial
   * and rational spline with respect to the outer function's control point
   * positions
   *
   * This function composes the Numerator only, so it can be reused to work with
   * rational-rational-spline-compositions
   */
  template <typename SplineType>
  constexpr auto ComposeNumeratorSensitivity(
      const SplineType &inner_function) const {
    return ComposeNumeratorHelper<true>(inner_function);
  }

  /**
   * @brief Compose the Numerator Function of a polynomial and rational spline
   *
   * This function composes the Numerator only, so it can be reused to work with
   * rational-rational-spline-compositions
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension_inner_spline,
                         decltype(ScalarType{} * ScalarRHS{}),
                         decltype(ScalarType{} * ScalarRHS{})>
  ComposeDenominatorSpline(
      const RationalBezierSpline<parametric_dimension_inner_spline,
                                 PointTypeRHS, ScalarRHS> &inner_function)
      const;

 public:
  /// Make Number Of Control Points available
  const IndexingType &GetNumberOfControlPoints() const {
    return number_of_control_points;
  };
  /// Make ScalarType publicly available
  using ScalarType_ = ScalarType;
  using IndexingType_ = IndexingType;
  using PhysicalPointType_ = PhysicalPointType;

  /// Offsets in Row based control point storage
  std::array<IndexingType, parametric_dimension> index_offsets{};
  /// List of all control points in "Row-based" order
  std::vector<PhysicalPointType_> control_points{};

  /// Make Parametric dimension publicly available
  static constexpr IndexingType kParametricDimensions = parametric_dimension;

  /*
   * Retrieve multidimensional indices
   *
   * @param local_indices single value index as control point index
   */
  constexpr std::array<IndexingType, parametric_dimension> LocalToGlobalIndex(
      const IndexingType &local_index) const;

  /*
   * Retrieve global scalar index
   *
   * @param global_idex Array index type as std::array
   */
  constexpr IndexingType GlobalToLocalIndex(
      const std::array<IndexingType, parametric_dimension> &global_index) const;

  /*
   * Calculate the new control point that result from the multiplication
   * between two bezier splines
   * Note : Using the notation in Gershons diss (eq. 2.13)
   */
  template <typename PhysicalPointLHS, typename ScalarLHS,
            typename PhysicalPointRHS, typename ScalarRHS, typename... T>
  constexpr void CombineControlPointsForProduct_(
      const BezierSpline<parametric_dimension, PhysicalPointLHS, ScalarLHS>
          &P_spline,
      const BezierSpline<parametric_dimension, PhysicalPointRHS, ScalarRHS>
          &Q_spline,
      const std::array<IndexingType, parametric_dimension> &ctpsIndex,
      const ScalarType factor, const T &...indices);

  /// Copy constructor
  constexpr BezierSpline(const BezierSpline &bezier_spline) = default;

  /// Empty constructor
  constexpr BezierSpline() = default;

  /// Empty constructor with degrees
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : degrees{deg} {
    number_of_control_points = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      number_of_control_points *= degrees[i] + 1;
    control_points.resize(number_of_control_points);
    UpdateIndexOffsets_();
  };

  /// Constructor with control point list
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> &deg,
      const std::vector<PhysicalPointType> &points)
      : degrees{deg}, control_points{points} {
    number_of_control_points = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      number_of_control_points *= degrees[i] + 1;

    UpdateIndexOffsets_();
    assert(number_of_control_points == points.size());
  };

  /// Move operator
  constexpr BezierSpline &operator=(BezierSpline &&rhs) = default;

  /// Move operator
  constexpr BezierSpline &operator=(const BezierSpline &rhs) = default;

  /// Getter for Degrees
  constexpr const std::array<std::size_t, parametric_dimension> &GetDegrees()
      const {
    return degrees;
  }

  /// Getter for Degrees
  constexpr std::array<std::size_t, parametric_dimension> GetDegrees() {
    return degrees;
  }

  /// Set Degrees
  constexpr void UpdateDegrees(
      const std::array<std::size_t, parametric_dimension> &new_degrees) {
    degrees = new_degrees;
    number_of_control_points = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      number_of_control_points *= degrees[i] + 1;
    control_points.resize(number_of_control_points);
    UpdateIndexOffsets_();
  }

  /// Access Control Point Vector directly (for conformity to rationals)
  constexpr const std::vector<PhysicalPointType> &GetControlPoints() const {
    return control_points;
  }

  /// Access Control Point Vector directly (for conformity to rationals)
  constexpr std::vector<PhysicalPointType> &GetControlPoints() {
    return control_points;
  }

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PhysicalPointType_ &ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PhysicalPointType_ &ControlPoint(const T... index);

  /// Retrieve single control point from local indices (as array)
  constexpr const PhysicalPointType_ &ControlPoint(
      const std::array<IndexingType, parametric_dimension> &index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr PhysicalPointType_ &ControlPoint(
      const std::array<IndexingType, parametric_dimension> &index);

  /// Order elevation along a specific parametric dimension
  constexpr BezierSpline &OrderElevateAlongParametricDimension(
      const IndexingType par_dim);

  /**
   * Determine the close form representation of the derivative
   *
   * The orders are set per parameteric dimension, meaning that {1, 2} will
   * result in the first order derivative along the first parametric axis
   * and the second order derivative alon the second parametric axis.
   *
   * The control points are computed as follows
   * \f[
   * \hat{C}_i^k=\frac{p!}{(p-k)!}\sum_{j=0}^k (-1)^{j+1}\binom{k}{j}C_{i+j}
   * \f]
   * where \f$\hat{C}_i^k\f$ denotes the control point of the close form of
   * spline the derivative, \f$ p \f$ is the degree of the original spline,
   * \f$k\f$ is the order of the derivative (in a given parametric direction)
   * \f$C_i\f$ are the original control points
   *
   * @param orders array denoting the orders along the respective parametric
   * axis
   */
  constexpr BezierSpline DerivativeWRTParametricDimension(
      const std::array<IndexingType, parametric_dimension> &orders) const;

  /// Evaluate the spline using the de Casteljau algorithm
  template <typename... T>
  constexpr PhysicalPointType_ Evaluate(const T &...par_coords) const {
    return Evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PhysicalPointType_ Evaluate(
      const PointTypeParametric_ &par_coords) const;

  /// Evaluate Basis Functions
  template <typename... T>
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctionContributions(const T &...par_coords) const {
    return BasisFunctionContributions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctionContributions(const PointTypeParametric_ &par_coords) const;

  /// Evaluate Basis Functions Unraveled using Cartesian Product of p-dim
  constexpr std::vector<ScalarType> BasisFunctions(
      const PointTypeParametric_ &par_coords) const;

  /// Evaluate Basis Functions Unraveled using Cartesian Product of p-dim
  template <typename... T>
  constexpr std::vector<ScalarType> BasisFunctions(
      const T &...par_coords) const {
    return BasisFunctions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions Derivatives
  constexpr std::vector<ScalarType> BasisFunctionsDerivatives(
      const PointTypeParametric_ &par_coords,
      const std::array<std::size_t, parametric_dimension> &nth_derivs) const;

  /// Evaluate Basis Functions Derivatives
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctionsDerivativeContributions(
      const PointTypeParametric_ &par_coords,
      const std::array<std::size_t, parametric_dimension> &nth_derivs) const;

  /// Evaluate the spline via the explicit precomputation of bernstein
  /// values
  constexpr PhysicalPointType_ ForwardEvaluate(
      const PointTypeParametric_ &par_coords) const;

  /// Evaluate the derivatives of a spline via the explicit precomputation of
  /// bernstein polynomial derivativs
  constexpr PhysicalPointType_ EvaluateDerivative(
      const PointTypeParametric_ &par_coords,
      const std::array<std::size_t, parametric_dimension> &nth_derivs) const;

  /// Evaluate the spline using explicit precomputation of bernstein values
  template <typename... T>
  constexpr PhysicalPointType_ ForwardEvaluate(const T &...par_coords) const {
    return ForwardEvaluate(PointTypeParametric_{par_coords...});
  }

  /// Addition of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator+(
      const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS> &rhs)
      const;

  /// Add two splines of same type
  constexpr BezierSpline &operator+=(BezierSpline rhs);

  /// Substraction of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator-(
      BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS> rhs) const;

  /// Add two splines of same type
  constexpr BezierSpline &operator-=(BezierSpline rhs);

  /// Check if two splines are equivalent
  constexpr bool operator==(const BezierSpline &rhs) const;

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator*(
      const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS> &rhs)
      const;

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
  ExtractDimension(const IndexingType &dimension) const;

  constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
  RaisePower(const IndexingType power) const;

  /// Multiplication with scalar
  constexpr BezierSpline &operator*=(const ScalarType &scalar);

  /// Multiplication with scalar RAII
  constexpr BezierSpline operator*(const ScalarType &scalar) const;

  /// Friend injection for reversed order
  friend constexpr BezierSpline operator*(const ScalarType &scalar,
                                          const BezierSpline &b) {
    return b * scalar;
  }

  /// Check if can be used for composition
  constexpr bool FitsIntoUnitCube() const;

  /// Addition
  constexpr BezierSpline &operator+=(const PhysicalPointType &point_shift);

  /// Addition
  constexpr BezierSpline operator+(const PhysicalPointType &point_shift) const;

  /// Friend injection Addition
  friend constexpr BezierSpline operator+(const PhysicalPointType &point_shift,
                                          const BezierSpline &original_spline) {
    return original_spline + point_shift;
  }

  /// Substraction
  constexpr BezierSpline operator-(const PhysicalPointType &point_shift) const;

  /// Substraction
  constexpr BezierSpline &operator-=(const PhysicalPointType &point_shift);

  /// Inversion
  constexpr BezierSpline operator-() const;

  /// Get maximum restricting corner of spline
  constexpr PhysicalPointType MaximumCorner() const;

  /// Get minimum restricting corner of spline
  constexpr PhysicalPointType MinimumCorner() const;

  /*
   * Reposition spline and scale dimensionwise
   *
   * Scales a spline along its different dimensions and transposes its position.
   * This can be used as an internal function when a spline is mapped into the
   * unit cube, but is also handy when the same operation is performed on a
   * group of splines where the scaling parameters and transposition vectors
   * stay constant for the entire group.
   *
   * @param transposition PointType describing the first corner of the nd cuboid
   * @param stretch PointType describing the second corner of the nd cuboid
   */
  constexpr BezierSpline &TransposeAndScale(
      const PhysicalPointType &transposition,
      const PhysicalPointType &scale_vector);

  /*
   * Fit into unit cube
   *
   * Takes spline and fits it into the unit cuboid (important for spline
   * composition)
   */
  constexpr BezierSpline &FitIntoUnitCube();

  /// Friend injection Substraction
  friend constexpr BezierSpline operator-(const PhysicalPointType &point_shift,
                                          const BezierSpline &original_spline) {
    return -(original_spline - point_shift);
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
  constexpr auto Compose(
      const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                         ScalarRHS> &inner_function) const;

  /*
   * Sensitivity of Functional Composition between two splines with respect to
   * outer spline's control-point position
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and the
   * function argument as the inner function. The result represents the
   * derivative of the functional composition with respect to the outer
   * geometries control point position, in the form of another Bezier Spline.
   * This works so long as the parametric dimension of the outer function
   * matches the physical dimension of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr auto ComposeSensitivity(
      const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                         ScalarRHS> &inner_function) const;

  /*
   * Functional Composition between a polynomial and rational spline
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and the
   * function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr auto Compose(
      const RationalBezierSpline<parametric_dimension_inner_spline,
                                 PointTypeRHS, ScalarRHS> &inner_function)
      const;
  /*
   * Sensitivity of Functional Composition between two splines with respect to
   * outer spline's control-point position
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and the
   * function argument as the inner function. The result represents the
   * derivative of the functional composition with respect to the outer
   * geometries control point position, in the form of another Bezier Spline.
   * This works so long as the parametric dimension of the outer function
   * matches the physical dimension of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr auto ComposeSensitivity(
      const RationalBezierSpline<parametric_dimension_inner_spline,
                                 PointTypeRHS, ScalarRHS> &inner_function)
      const;

  /*
   * Composition between mutliple splines from a spline group
   *
   * Performes a composition between multple splines, which can be used to
   * construct microstructures. After the return group is instantiated, the
   * composition is performed elementwise.
   */
  template <typename SplineType>
  constexpr auto Compose(
      const std::vector<SplineType> &inner_function_group) const;

  /*
   * Split the Bezier Spline into two distinct subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over two splines.
   */
  constexpr std::vector<BezierSpline> SplitAtPosition(
      const ScalarType &splitting_plane,
      const IndexingType splitting_dimension = 0) const;

  /*
   * Split the Bezier Spline into several subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over several splines.
   */
  constexpr std::vector<BezierSpline> SplitAtPosition(
      const std::vector<ScalarType> &splitting_planes,
      const IndexingType splitting_dimension = 0) const;
};

#include "bezman/src/bezier_spline.inc"

}  // namespace bezman

#endif  // SRC_BEZIER_SPLINE_HPP
