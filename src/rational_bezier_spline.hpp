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

#ifndef SRC_RATIONAL_BEZIER_SPLINE_HPP
#define SRC_RATIONAL_BEZIER_SPLINE_HPP

#include <array>
#include <vector>

#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"

namespace bezman {

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class BezierSpline;

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
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;
  using PolynomialBezier =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;

#ifndef BEZMAN_DYNAMIC
  using PolynomialScalarBezier =
      BezierSpline<parametric_dimension, ScalarType, ScalarType>;
#else
  using PolynomialScalarBezier = PolynomialBezier;
#endif
  /*!
   * Numerator Spline, with control points multiplied with their
   * respective weights
   */
  PolynomialBezier weighted_spline_;

  /*!
   * Denominator Spline contains the weight function
   */
  PolynomialScalarBezier weight_function_;

  /**
   * @brief  Wrapper function that checks if the Numerator and Denominator
   * splines are still compatible after operations are performed
   */
  constexpr bool CheckSplineCompatibility() const {
    bool compatible =
        (weighted_spline_.GetDegrees() == weight_function_.GetDegrees());
    return compatible;
  }

  // Friend declarations
  template <std::size_t parent_parametric_dimension,
            typename ParentPhysicalPointType, typename ParentScalarType>
  friend class BezierSpline;
  template <std::size_t parent_parametric_dimension,
            typename ParentPhysicalPointType, typename ParentScalarType>
  friend class RationalBezierSpline;

 public:
  /// Make ScalarType publicly available
  using ScalarType_ = ScalarType;
  using PhysicalPointType_ = PhysicalPointType;
#ifdef BEZMAN_DYNAMIC
  using ControlPointsType_ = typename PolynomialBezier::ControlPointsType_;
  using ReturnVectorType_ = typename PolynomialBezier::ReturnVectorType_;
#endif

  /// Make Parametric dimension publicly available
  static constexpr IndexingType kParametricDimensions = parametric_dimension;

  /// Copy constructor
  constexpr RationalBezierSpline(const RationalBezierSpline& bezier_spline) =
      default;

  /// Transform Bezier into rational Bezier
  constexpr RationalBezierSpline(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
          bezier_spline)
      : weighted_spline_{bezier_spline} {
#ifndef BEZMAN_DYNAMIC
    // Initialize new weight function spline
    weight_function_ =
        BezierSpline<parametric_dimension, ScalarType, ScalarType>{
            // Degrees can be taken from weighted spline
            weighted_spline_.GetDegrees(),
            // Generate control weights (just 1 for every point)
            std::vector<ScalarType>(
                // Size of control-point vector
                weighted_spline_.GetNumberOfControlPoints(),
                // Value
                static_cast<ScalarType>(1.))};
#else
    weight_function_ =
        BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>{
            weighted_spline_.GetDegrees(),
            ControlPointsType_{std::vector<ScalarType>(
                                   weighted_spline_.GetNumberOfControlPoints(),
                                   static_cast<ScalarType>(1.)),
                               1}};
#endif
    assert(CheckSplineCompatibility());
  }

  /// Constructor with two spline - weight and weighted spline
  constexpr RationalBezierSpline(
      const PolynomialBezier& r_weighted_spline_,
      const PolynomialScalarBezier& r_weight_function_)
      : weighted_spline_{r_weighted_spline_},
        weight_function_{r_weight_function_} {
    assert(CheckSplineCompatibility());
  }

  /// Empty constructor
  constexpr RationalBezierSpline() = default;

#ifndef BEZMAN_DYNAMIC
  /// Constructor with weights
  constexpr RationalBezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : weighted_spline_{deg}, weight_function_{deg} {};
#else
  /// Constructor with weights
  constexpr RationalBezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const IndexingType dimension)
      : weighted_spline_{deg, dimension}, weight_function_{deg, dimension} {};
#endif

  /// Constructor with control point list
  constexpr RationalBezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const ControlPointsType_ control_point_vector,
      const std::vector<ScalarType> weights)
#ifndef BEZMAN_DYNAMIC
      : weighted_spline_{deg}, weight_function_{deg, weights} {
    assert(control_point_vector.size() == weights.size());
    assert(weighted_spline_.GetNumberOfControlPoints() == weights.size());
    for (std::size_t i_point{0};
         i_point < weighted_spline_.control_points.size(); i_point++) {
      weighted_spline_.control_points[i_point] =
          control_point_vector[i_point] * weights[i_point];
    }
#else
  {
    const auto cp_vector_size = control_point_vector.size();
    const auto w_size = weights.size();
    const auto dim = control_point_vector.GetDimension();

    if (cp_vector_size != w_size) {
      throw std::runtime_error(
          "bezman::RationalBezierSpline - control points and weight size "
          "mismatch.");
    }
    ControlPointsType_ weighted_vertices{};
    auto& weighted_vertices_vector = weighted_vertices.GetVertices();
    weighted_vertices_vector.reserve(cp_vector_size);

    auto cp_vector_iter = control_point_vector.GetVertices().cbegin();
    for (std::size_t i_point{0}; i_point < cp_vector_size; ++i_point) {
      const auto& w = weights[i_point];
      for (std::size_t i_dim{}; i_dim < dim; ++i_dim) {
        weighted_vertices_vector.emplace_back(*(cp_vector_iter++) * w);
      }
    }

    weighted_vertices.SetDimension(dim);

    weighted_spline_ = PolynomialBezier{deg, weighted_vertices};
    weight_function_ = PolynomialScalarBezier{deg, Vertices{weights, 1}};
#endif
    assert(CheckSplineCompatibility());
  }

  /// Move operator
  constexpr RationalBezierSpline& operator=(RationalBezierSpline&& rhs) =
      default;

  /// Move operator
  constexpr RationalBezierSpline& operator=(const RationalBezierSpline& rhs) =
      default;

  /// Getter for Degrees
  constexpr const std::array<std::size_t, parametric_dimension>& GetDegrees()
      const {
    return weighted_spline_.GetDegrees();
  }

  /// Getter for Number Of Control Points
  const IndexingType& GetNumberOfControlPoints() const {
    return weighted_spline_.GetNumberOfControlPoints();
  };

  /// Access Control Point Vector directly
  constexpr const ControlPointsType_& GetWeightedControlPoints() const {
    return weighted_spline_.control_points;
  }

  /// Access Control Point Vector directly
  constexpr ControlPointsType_& GetWeightedControlPoints() {
    return weighted_spline_.control_points;
  }

  /// Access Weights Vector directly
  constexpr const std::vector<ScalarType>& GetWeights() const {
    return weight_function_.control_points.GetVertices();
  }

  /// Access Weights Vector directly
  constexpr std::vector<ScalarType>& GetWeights() {
    return weight_function_.control_points.GetVertices();
  }

  /// Set Degrees
  constexpr void UpdateDegrees(
      const std::array<std::size_t, parametric_dimension>& new_degrees) {
    weight_function_.UpdateDegrees(new_degrees);
    weighted_spline_.UpdateDegrees(new_degrees);
    assert(CheckSplineCompatibility());
    return;
  }

  // @NOTE These functions all pass arguments and perform very limited
  // operations

  // Get the weight function spline as const reference
  constexpr const PolynomialScalarBezier& GetWeightFunctionSpline() const {
    return weight_function_;
  }

#ifndef BEZMAN_DYNAMIC
  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PhysicalPointType_ ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr PhysicalPointType_ ControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted eighted control point from local indices
  template <typename... T>
  constexpr const PhysicalPointType_& WeightedControlPoint(
      const T... index) const;

  /// Retrieve single weighted control point from local indices
  template <typename... T>
  constexpr PhysicalPointType_& WeightedControlPoint(const T... index);

  /// Retrieve single weighted control point from local indices (as array)
  constexpr const PhysicalPointType_& WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted control point from local indices (as array)
  constexpr PhysicalPointType_& WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index);

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr const ScalarType& Weight(const T... index) const;

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr ScalarType& Weight(const T... index);

  /// Retrieve single weight from local indices (as array)
  constexpr const ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weight from local indices (as array)
  constexpr ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index);
#else
  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr auto ControlPoint(const T... index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr auto ControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted eighted control point from local indices
  template <typename... T>
  constexpr const auto WeightedControlPoint(const T... index) const;

  /// Retrieve single weighted control point from local indices
  template <typename... T>
  constexpr auto WeightedControlPoint(const T... index);

  /// Retrieve single weighted control point from local indices (as array)
  constexpr const auto WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weighted control point from local indices (as array)
  constexpr auto WeightedControlPoint(
      const std::array<IndexingType, parametric_dimension>& index);

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr const ScalarType& Weight(const T... index) const;

  /// Retrieve single weight from local indices
  template <typename... T>
  constexpr ScalarType& Weight(const T... index);

  /// Retrieve single weight from local indices (as array)
  constexpr const ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single weight from local indices (as array)
  constexpr ScalarType& Weight(
      const std::array<IndexingType, parametric_dimension>& index);
#endif
  //------------------- Essential Operations

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr auto Evaluate(const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using the de Casteljau algorithm
  template <typename... T>
  constexpr auto Evaluate(const T&... par_coords) const {
    return Evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions in the form (w_i N_i / (sum(w_j N_j)))
  template <typename... T>
  constexpr std::vector<ScalarType> BasisFunctions(
      const T&... par_coords) const {
    return BasisFunctions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions
  constexpr std::vector<ScalarType> BasisFunctions(
      const PointTypeParametric_& par_coords) const;

  /**
   * @brief Evaluate non weighted basis functions
   *
   * returns basis function values in the form (N_i / (sum(w_j N_j)))
   * They are required in the use of weighted control points to avoid
   * multiplication followed by division of the same coefficients. For internal
   * use only
   */
  template <typename... T>
  constexpr std::vector<ScalarType> UnweightedBasisFunctions(
      const T&... par_coords) const {
    return UnweightedBasisFunctions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions
  constexpr std::vector<ScalarType> UnweightedBasisFunctions(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate Basis Functions
  template <typename... T>
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctionContributions(const T&... par_coords) const {
    return BasisFunctionContributions(PointTypeParametric_{par_coords...});
  }

  /// Evaluate Basis Functions without respecting weights
  constexpr std::array<std::vector<ScalarType>, parametric_dimension>
  BasisFunctionContributions(const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline via the explicit precomputation of bernstein values
  constexpr PhysicalPointType_ ForwardEvaluate(
      const PointTypeParametric_& par_coords) const;

  /// Evaluate the spline using explicit precomputation of bernstein values
  template <typename... T>
  constexpr PhysicalPointType_ ForwardEvaluate(const T&... par_coords) const {
    return ForwardEvaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate the derivatives of a spline using Leibnitz' rule
  constexpr auto EvaluateDerivative(
      const PointTypeParametric_& par_coords,
      const std::array<std::size_t, parametric_dimension>& nth_derivs) const;

  /// Evaluate the derivatives of a spline using Leibnitz' rule
  constexpr std::vector<ScalarType> BasisFunctionsDerivatives(
      const PointTypeParametric_& par_coords,
      const std::array<std::size_t, parametric_dimension>& nth_derivs) const;

  //-------------------  Refinement Strategies

  /// Order elevation along a specific parametric dimension
  constexpr RationalBezierSpline& OrderElevateAlongParametricDimension(
      const IndexingType par_dim);

  /*
   * Split the Bezier Spline into two distinct subdivisions
   *
   * Splits the Spline along a specific dimension and returns a group
   * representing the same domain over two splines.
   */
  constexpr BezierGroup<RationalBezierSpline> SplitAtPosition(
      const ScalarType& splitting_plane,
      const IndexingType splitting_dimension = 0) const;

  /*
   * Split the Bezier Spline into distinct subdivisions along vector
   *
   * Splits the Spline along a specific dimension at entries specified within a
   * vector and returns a group representing the same domain over n+1 splines
   */
  constexpr BezierGroup<RationalBezierSpline> SplitAtPosition(
      const std::vector<ScalarType>& splitting_planes,
      const IndexingType splitting_dimension = 0) const;

  /**
   * @brief  Derivative along a specific parametric dimension
   *
   * This function is implemented using the Quotient-rule for derivatives
   */
  constexpr RationalBezierSpline DerivativeWRTParametricDimension(
      const IndexingType par_dim) const;

  //-------------------  Operator overloads

  /// Multiplication with scalar
  constexpr RationalBezierSpline& operator*=(const ScalarType& scalar);

  /// Multiplication with scalar RAII
  constexpr RationalBezierSpline operator*(const ScalarType& scalar) const;

  /// Friend injection for reversed order
  friend constexpr RationalBezierSpline operator*(
      const ScalarType& scalar, const RationalBezierSpline& b) {
    return b * scalar;
  }

  /// Addition
  constexpr RationalBezierSpline& operator+=(
      const PhysicalPointType& point_shift);

  /// Addition
  constexpr RationalBezierSpline operator+(
      const PhysicalPointType& point_shift) const;

  /// Friend injection Addition
  friend constexpr RationalBezierSpline operator+(
      const PhysicalPointType& point_shift,
      const RationalBezierSpline& original_spline) {
    return original_spline + point_shift;
  }

  /// Substraction
  constexpr RationalBezierSpline operator-(
      const PhysicalPointType& point_shift) const;

  /// Substraction
  constexpr RationalBezierSpline& operator-=(
      const PhysicalPointType& point_shift);

  /// Inversion
  constexpr RationalBezierSpline operator-() const;

  /// Addition of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator+(
      const RationalBezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
          rhs) const;

  /// Add two splines of same type
  constexpr RationalBezierSpline& operator+=(const RationalBezierSpline& rhs);

  /// Substraction of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator-(
      const RationalBezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
          rhs) const;

  /// Subtract two splines of same type
  constexpr RationalBezierSpline& operator-=(const RationalBezierSpline& rhs);

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr auto operator*(
      const RationalBezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
          rhs) const;

  //------------------- IRIT-related operations

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
  ExtractDimension(const std::size_t& dimension) const;

  /*!
   * Functional Composition between two splines
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and
   * the function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr auto Compose(
      const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                         ScalarRHS>& inner_function) const {
    return RationalBezierSpline<parametric_dimension_inner_spline,
                                PhysicalPointType,
                                decltype(ScalarType_{} * ScalarRHS{})>{
        // Numerator composition
        weighted_spline_.Compose(inner_function),
        // Denominator composition
        weight_function_.Compose(inner_function)};
  }

  /*!
   * Sensitivity to Functional Composition between two splines
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
                         ScalarRHS>& inner_function) const;

  /*!
   * Sensitivity to Functional Composition between two splines
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
                                 PointTypeRHS, ScalarRHS>& inner_function)
      const;

  /*!
   * Functional Composition between two rational
   *
   * Compose two splines, taking the (*this) spline as the outer funtion and
   * the function argument as the inner function. This works so long as the
   * parametric dimension of the outer function matches the physical dimension
   * of the inner function.
   *
   * For two rational splines, the weight function cancels out, see
   * BezierSpline::ComposeDenominator
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr RationalBezierSpline<parametric_dimension_inner_spline,
                                 PhysicalPointType,
                                 decltype(ScalarType_{} * ScalarRHS{})>
  Compose(const RationalBezierSpline<parametric_dimension_inner_spline,
                                     PointTypeRHS, ScalarRHS>& inner_function)
      const {
    return RationalBezierSpline<parametric_dimension_inner_spline,
                                PhysicalPointType,
                                decltype(ScalarType_{} * ScalarRHS{})>{
        // Numerator composition
        weighted_spline_.ComposeNumeratorSpline(inner_function),
        // Denominator composition
        weight_function_.ComposeNumeratorSpline(inner_function)};
  }

  /// Check if can be used for composition
  constexpr bool FitsIntoUnitCube() const {
    return weighted_spline_.FitsIntoUnitCube();
  };

  /// Get maximum restricting corner of spline
  constexpr auto MaximumCorner() const {
    return weighted_spline_.MaximumCorner();
  };

  /// Get minimum restricting corner of spline
  constexpr auto MinimumCorner() const {
    return weighted_spline_.MinimumCorner();
  };
};

#include "bezman/src/rational_bezier_spline.inc"

}  // namespace bezman

#endif  // SRC_RATIONAL_BEZIER_SPLINE_HPP
