#ifndef SRC_BEZIER_HPP
#define SRC_BEZIER_HPP

#include <array>
#include <cassert>
#include <numeric>
#include <vector>

#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/fastbinomialcoefficient.hpp"

namespace beziermanipulation {

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
 private:
  using IndexingType = std::size_t;
  using PointTypePhysical_ = PhysicalPointType;
  using PointTypeParametric_ = Point<parametric_dimension, ScalarType>;

  constexpr void update_index_offsets();

 public:
  using ScalarType_ = ScalarType;

  /// Polynomial degrees
  std::array<IndexingType, parametric_dimension> degrees;
  /// Offsets in Row based control point storage
  std::array<IndexingType, parametric_dimension> index_offsets;
  /// Number of control points
  IndexingType NumberOfControlPoints{};
  /// List of all control points in "Row-based" order
  std::vector<PointTypePhysical_> control_points{};

  /*
   * Retrieve individual indices
   *
   * @param local_indices single value index as control point index
   */
  constexpr std::array<IndexingType, parametric_dimension>
  local_to_global_index(const IndexingType& local_index);

  /*
   * Calculate the new control point that result from the multiplication
   * between two bezier splines
   * Note : Using the notation in Gershons diss (eq. 2.13)
   */
  template <typename PhysicalPointLHS, typename ScalarLHS,
            typename PhysicalPointRHS, typename ScalarRHS, typename... T>
  constexpr void product_combine_control_points(
      const BezierSpline<parametric_dimension, PhysicalPointLHS, ScalarLHS>&
          P_spline,
      const BezierSpline<parametric_dimension, PhysicalPointRHS, ScalarRHS>&
          Q_spline,
      const std::array<IndexingType, parametric_dimension>& ctpsIndex,
      const ScalarType factor, const T&... indices);

  /// Copy constructor
  constexpr BezierSpline(const BezierSpline& bezier_spline) = default;

  /// Empty constructor
  constexpr BezierSpline() = default;

  /// Empty constructor
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : degrees{deg} {
    NumberOfControlPoints = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      NumberOfControlPoints *= degrees[i] + 1;
    control_points.resize(NumberOfControlPoints);
    update_index_offsets();
  };

  /// Constructor with control point list
  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const std::vector<PointTypePhysical_> points)
      : degrees{deg}, control_points{points} {
    NumberOfControlPoints = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      NumberOfControlPoints *= degrees[i] + 1;

    update_index_offsets();
    assert(NumberOfControlPoints == points.size());
  };

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& control_point(const T... index) const {
    static_assert(sizeof...(T) == parametric_dimension,
                  "Unspecified number of indices.");
    unsigned int c_i{0}, i{};
    ((c_i += index_offsets[i++] * index), ...);
    return control_points[c_i];
  }

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& control_point(const T... index) {
    static_assert(sizeof...(T) == parametric_dimension,
                  "Unspecified number of indices.");
    unsigned int c_i{0}, i{};
    ((c_i += index_offsets[i++] * index), ...);
    return control_points[c_i];
  }

  /// Retrieve single control point from local indices (as array)
  constexpr const PointTypePhysical_& control_point(
      const std::array<IndexingType, parametric_dimension>& index) const {
    unsigned int c_i{0}, i{};
    for (unsigned int i{}; i < parametric_dimension; i++) {
      c_i += index_offsets[i] * index[i];
    }
    return control_points[c_i];
  }

  /// Retrieve single control point from local indices (as array)
  constexpr PointTypePhysical_& control_point(
      const std::array<IndexingType, parametric_dimension>& index) {
    unsigned int c_i{0}, i{};
    for (unsigned int i{}; i < parametric_dimension; i++) {
      c_i += index_offsets[i] * index[i];
    }
    return control_points[c_i];
  }

  /// Order elevation along a specific parametric dimension
  constexpr BezierSpline& order_elevate_along_parametric_dimension(
      const IndexingType par_dim);

  template <typename... T>
  constexpr PointTypePhysical_ evaluate(const T&... par_coords) const {
    return (*this).evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ evaluate(
      const PointTypeParametric_& par_coords) const;

  /// Addition of Two Splines resulting in a new spline that describes the
  /// pointwise addition of the two Beziers
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} + PointTypeRHS{}),
                         decltype(ScalarType_{} + ScalarRHS{})>
  operator+(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const;

  /// Check if two splines are equivalent
  constexpr bool operator==(const BezierSpline& rhs) const {
    // Check if degrees fit and if control_points are the same
    if (rhs.degrees == degrees) {
      return rhs.control_points == control_points;
    } else {
      return false;
    }
  }

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} * PointTypeRHS{}),
                         decltype(ScalarType_{} * ScalarRHS{})>
  operator*(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const;

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
  ExtractDimension(unsigned int dimension) const {
    assert(dimension < PointTypePhysical_::kSpatialDimension);
    BezierSpline<parametric_dimension, ScalarType, ScalarType> extracted_spline(
        degrees);
    for (std::size_t i{}; i < NumberOfControlPoints; i++) {
      extracted_spline.control_points[i] = control_points[i][dimension];
    }
    return extracted_spline;
  }

  constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
  power(const unsigned int power) const {
    // @TODO use log2(power) algorithm that squares the result to minimize
    // multiplications
    static_assert(std::is_scalar_v<PhysicalPointType>,
                  "Only Scalar-type Splines can be raised to a power.");
    assert(
        ("Not implemented, as raising to 0 would be inefficient.", power > 0));
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
        power_spline{(*this)};
    for (int i{1}; i < power; i++) {
      power_spline = power_spline * (*this);
    }
    return power_spline;
  }

  /// Multiplication with scalar
  constexpr BezierSpline operator*(const ScalarType& scalar) const {
    BezierSpline scaled_spline{(*this)};
    for (std::size_t i{}; i < NumberOfControlPoints; i++) {
      scaled_spline.control_points[i] *= scalar;
    }
  };

  /// Frind injection
  friend constexpr BezierSpline operator*(const ScalarType& scalar,
                                          const BezierSpline& b) {
    return b * scalar;
  }

  /// Compose two splines
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                         decltype(ScalarType_{} * ScalarRHS{})>
  compose(const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                             ScalarRHS>& inner_function) const {
    /// Start the composition with the current spline
    // Initialize return value
    using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});
    const IndexingType sum_of_degrees_outer_spline =
        std::accumulate(degrees.begin(), degrees.end(), 0);
    // New degrees
    std::array<IndexingType, parametric_dimension_inner_spline> new_degrees{
        inner_function.degrees};
    for (IndexingType i{}; i < parametric_dimension_inner_spline; i++)
      new_degrees[i] *= sum_of_degrees_outer_spline;

    BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                 ScalarReturnT>
        composition{new_degrees};

    // Check Dimensions
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");

    // Loop over the parametric dimensions of the outer function
    for (std::size_t i_outer_parametric_dimension{};
         i_outer_parametric_dimension < parametric_dimension;
         i_outer_parametric_dimension++) {
      // Extract the current dimension of the spline
      const auto inner_spline_xi =
          inner_function.ExtractDimension(i_outer_parametric_dimension);
      const IndexingType n_outer_degree = degrees[i_outer_parametric_dimension];

      // First -> Calculate the individual Basis functions, which are now
      // represented as splines
      for (std::size_t i_basis_outer_spline{};
           i_basis_outer_spline <= n_outer_degree; i_basis_outer_spline++) {
        // Composed Basis function
        std::vector<ScalarReturnT> scaling_factors{};
        if ((i_basis_outer_spline == 0) ||
            (i_basis_outer_spline == n_outer_degree)) {
          // Boundary cases
          const auto B_i = inner_spline_xi.power(n_outer_degree);
          scaling_factors = B_i.control_points;
        } else {
          const auto B_i =
              utils::FastBinomialCoefficient::choose(n_outer_degree,
                                                     i_basis_outer_spline) *
              inner_spline_xi.power(i_basis_outer_spline) *
              inner_spline_xi.power(n_outer_degree - i_basis_outer_spline);
          scaling_factors = B_i.control_points;
        }
        // for (IndexingType i_par_dim{}; i_par_dim < parametric_dimension) }
      }
    }
    return composition;
  }
};

#include "bezierManipulation/src/bezier_spline.inc"

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_HPP
