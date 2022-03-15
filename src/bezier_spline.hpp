#ifndef SRC_BEZIER_SPLINE_HPP
#define SRC_BEZIER_SPLINE_HPP

#include <array>
#include <cassert>
#include <numeric>
#include <vector>

#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/fastbinomialcoefficient.hpp"

namespace beziermanipulation {

// Forward declaration for later use
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
class BezierSplineGroup;

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
  local_to_global_index(const IndexingType& local_index) const;

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

  /// Empty constructor with degrees
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

  /// Move operator
  constexpr BezierSpline& operator=(BezierSpline&& rhs) = default;

  /// Move operator
  constexpr BezierSpline& operator=(const BezierSpline& rhs) = default;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr const PointTypePhysical_& control_point(const T... index) const;

  /// Retrieve single control point from local indices
  template <typename... T>
  constexpr PointTypePhysical_& control_point(const T... index);

  /// Retrieve single control point from local indices (as array)
  constexpr const PointTypePhysical_& control_point(
      const std::array<IndexingType, parametric_dimension>& index) const;

  /// Retrieve single control point from local indices (as array)
  constexpr PointTypePhysical_& control_point(
      const std::array<IndexingType, parametric_dimension>& index);

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
  constexpr bool operator==(const BezierSpline& rhs) const;

  /// Multiplication of two splines similar to pointwise product
  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} * PointTypeRHS{}),
                         decltype(ScalarType_{} * ScalarRHS{})>
  operator*(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const;

  /// Extract single coordinate spline
  constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
  ExtractDimension(unsigned int dimension) const;

  constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
  power(const unsigned int power) const;

  /// Multiplication with scalar
  constexpr BezierSpline& operator*=(const ScalarType& scalar);

  /// Multiplication with scalar RAII
  constexpr BezierSpline operator*(const ScalarType& scalar) const;

  /// Friend injection for reversed order
  friend constexpr BezierSpline operator*(const ScalarType& scalar,
                                          const BezierSpline& b) {
    return b * scalar;
  }

  /// Check if can be used for composition
  constexpr bool fits_unit_cube() const;

  /// Addition
  constexpr BezierSpline& operator+=(const PhysicalPointType& point_shift);

  /// Addition
  constexpr BezierSpline operator+(const PhysicalPointType& point_shift) const;

  /// Friend injection Addition
  friend constexpr BezierSpline operator+(const PhysicalPointType& point_shift,
                                          const BezierSpline& original_spline) {
    return original_spline + point_shift;
  }

  /// Substraction
  constexpr BezierSpline operator-(const PhysicalPointType& point_shift) const;

  /// Substraction
  constexpr BezierSpline& operator-=(const PhysicalPointType& point_shift);

  /// Inversion
  constexpr BezierSpline operator-() const;

  /// Get maximum restricting corner of spline
  constexpr PhysicalPointType maximum() const;

  /// Get minimum restricting corner of spline
  constexpr PhysicalPointType minimum() const;

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
  constexpr BezierSpline& transpose_and_scale(
      const PhysicalPointType& transposition,
      const PhysicalPointType& scale_vector);

  /*
   * Fit into unit cube
   *
   * Takes spline and fits it into the unit cuboid (important for spline
   * composition)
   */
  constexpr BezierSpline& fit_to_unit_cube();

  /// Friend injection Substraction
  friend constexpr BezierSpline operator-(const PhysicalPointType& point_shift,
                                          const BezierSpline& original_spline) {
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
  constexpr BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                         decltype(ScalarType_{} * ScalarRHS{})>
  compose(const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                             ScalarRHS>& inner_function) const;

  /*
   * Composition between mutliple splines from a spline group
   *
   * Performes a composition between multple splines, which can be used to
   * construct microstructures. After the return group is instantiated, the
   * composition is performed elementwise.
   */
  template <std::size_t parametric_dimension_inner_spline,
            typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSplineGroup<parametric_dimension_inner_spline,
                              PhysicalPointType,
                              decltype(ScalarType_{} * ScalarRHS{})>
  compose(
      const BezierSplineGroup<parametric_dimension_inner_spline, PointTypeRHS,
                              ScalarRHS>& inner_function_group) const;
  
};

#include "bezierManipulation/src/bezier_spline.inc"

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_SPLINE_HPP
