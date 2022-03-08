#ifndef SRC_BEZIER_HPP
#define SRC_BEZIER_HPP

#include <array>
#include <cassert>
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

  constexpr void update_index_offsets() {
    index_offsets[0] = 1;
    for (unsigned int i{1}; i < parametric_dimension; i++)
      index_offsets[i] = index_offsets[i - 1] * (degrees[i - 1] + 1);
  }

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
  local_to_global_index(const IndexingType& local_index) {
    std::array<IndexingType, parametric_dimension> indexList{};
    for (unsigned int i{0}; i < parametric_dimension; i++) {
      indexList[i] = (local_index / index_offsets[i]) % (degrees[i] + 1);
    }
    return indexList;
  }

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
      const ScalarType factor, const T&... indices) {
    // Some constant indices and degrees
    const int depth = sizeof...(indices);
    const int k = ctpsIndex[depth];
    const int m = P_spline.degrees[depth];
    const int n = Q_spline.degrees[depth];

    // Loop over current parametric domain
    for (int i{std::max(0, k - n)}; i <= std::min(k, m); i++) {
      // Calculate Factor
      const ScalarType lFactor =
          static_cast<ScalarType>(
              utils::FastBinomialCoefficient::choose(m, i) *
              utils::FastBinomialCoefficient::choose(n, k - i)) /
          static_cast<ScalarType>(
              utils::FastBinomialCoefficient::choose(m + n, k));

      // Now decide if continue recursion
      if constexpr ((depth + 1) == parametric_dimension) {
        const std::array<IndexingType, parametric_dimension> ind_lhs{
            static_cast<IndexingType>(indices)...,
            static_cast<IndexingType>(i)};
        std::array<IndexingType, parametric_dimension> ind_rhs{};
        for (unsigned int j{}; j < parametric_dimension; j++) {
          ind_rhs[j] = ctpsIndex[j] - ind_lhs[j];
        }
        (*this).control_point(ctpsIndex) += P_spline.control_point(ind_lhs) *
                                            Q_spline.control_point(ind_rhs) *
                                            factor * lFactor;
      } else {
        product_combine_control_points(P_spline, Q_spline, ctpsIndex,
                                       factor * lFactor, indices..., i);
      }
    }
  }

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
      const IndexingType par_dim) {
    // Calculate index Offsets to facilitate working on 1D array
    const unsigned int n_starting_points =
        (NumberOfControlPoints / (degrees[par_dim] + 1));
    const unsigned int starting_point_offset =
        index_offsets[par_dim] * (degrees[par_dim] + 1);
    const int starting_points_per_group = index_offsets[par_dim];
    const int n_groups = n_starting_points / starting_points_per_group;

    // Resize the CTPS vector accordingly
    NumberOfControlPoints =
        NumberOfControlPoints / (degrees[par_dim] + 1) * (degrees[par_dim] + 2);
    control_points.resize(NumberOfControlPoints);
    degrees[par_dim]++;

    // Local Counter
    unsigned int global_index = NumberOfControlPoints - 1;

    // Precalculations
    const ScalarType_ inverse_factor =
        static_cast<ScalarType_>(1) /
        static_cast<ScalarType_>(degrees[par_dim]);
    const IndexingType variable_offset_factor =
        index_offsets[par_dim] * (degrees[par_dim] - 1);

    // Vector is calculated from back to front, to hinder overwrite
    for (int group_index{n_groups - 1}; group_index >= 0; group_index--) {
      // Local variables
      const unsigned int first_index_in_group =
          group_index * (starting_point_offset);
      IndexingType i = degrees[par_dim] - 1;

      // Fix the last entry for element in the group
      for (int index_in_group{starting_points_per_group - 1};
           index_in_group >= 0; index_in_group--) {
        control_points[global_index] =
            control_points[first_index_in_group + index_in_group +
                           variable_offset_factor];
        global_index--;
      }

      // Interpolate for all but the first points in the vector, constantly
      // decreasing the counter index algorithm found in
      // https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/spline/Bezier/bezier-elev.html
      // The complex indexation is a result from the row based storage of the
      // control point positions
      for (IndexingType i{degrees[par_dim] - 1}; i > 0; i -= 1) {
        for (int index_in_group{starting_points_per_group - 1};
             index_in_group >= 0; index_in_group--) {
          const ScalarType_ factor =
              static_cast<ScalarType_>(i) * inverse_factor;

          control_points[global_index] =
              control_points[first_index_in_group + index_in_group +
                             index_offsets[par_dim] * i] *
                  (1 - factor) +
              control_points[first_index_in_group + index_in_group +
                             index_offsets[par_dim] * (i - 1)] *
                  factor;
          global_index--;
        }
      }

      // Fixate the first entry along each parametric dimension
      for (int index_in_group{starting_points_per_group - 1};
           index_in_group >= 0; index_in_group--) {
        // No we actually start the algorithm
        control_points[global_index] =
            control_points[first_index_in_group + index_in_group];
        global_index--;
      }
    }
    update_index_offsets();
    return (*this);
  }

  template <typename... T>
  constexpr PointTypePhysical_ evaluate(
      const T&... par_coords) const {
    return (*this).evaluate(PointTypeParametric_{par_coords...});
  }

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ evaluate(const PointTypeParametric_& par_coords)
      const {  // Work on copy of control_point
    std::vector<PointTypePhysical_> ctps_copy{control_points};
    IndexingType ctps_to_consider = NumberOfControlPoints;

    for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
      ScalarType_ factor = par_coords[par_dim];
      ScalarType_ inv_factor = 1. - par_coords[par_dim];

      ctps_to_consider /= degrees[par_dim] + 1;

      // For every starting position
      for (IndexingType start{0}; start < ctps_to_consider; start++) {
        const auto offset = index_offsets[par_dim + 1] * start;
        const auto step_width = index_offsets[par_dim];

        for (IndexingType i{0}; i <= degrees[par_dim]; i++) {
          for (IndexingType j{0}; j < degrees[par_dim] - i; j++) {
            ctps_copy[j * step_width + offset] =
                ctps_copy[j * step_width + offset] * inv_factor +
                ctps_copy[(j + 1) * step_width + offset] * factor;
          }
        }
      }
    }
    return ctps_copy[0];
  }

  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} + PointTypeRHS{}),
                         decltype(ScalarType_{} + ScalarRHS{})>
  operator+(const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>&
                rhs) const {
    // Initialize return value
    using PointTypeReturnT = decltype(PhysicalPointType{} + PointTypeRHS{});
    using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});

    BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
        return_spline(degrees, control_points);

    // Check if the right hand side requires a copy as it should not be
    // altered for this purpose
    bool rhs_needs_copy = false;
    for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
      rhs_needs_copy =
          rhs_needs_copy || (degrees[par_dim] > rhs.degrees[par_dim]);
    }

    // Increase the order of the copied spline to be greater or equal to the
    // RHSs order
    for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
      while (rhs.degrees[par_dim] > return_spline.degrees[par_dim]) {
        return_spline.order_elevate_along_parametric_dimension(par_dim);
      }
    }
    if (rhs_needs_copy) {
      // use commutativity of addition to create a copy of rhs
      return rhs + return_spline;
    } else {
      for (IndexingType i_ctps{}; i_ctps < return_spline.NumberOfControlPoints;
           i_ctps++) {
        return_spline.control_points[i_ctps] += rhs.control_points[i_ctps];
      }
      return return_spline;
    }
  }

  constexpr bool operator==(const BezierSpline& rhs) const {
    // Check if degrees fit and if control_points are the same
    if (rhs.degrees == degrees) {
      return rhs.control_points == control_points;
    } else {
      return false;
    }
  }

  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} * PointTypeRHS{}),
                         decltype(ScalarType_{} * ScalarRHS{})>
  operator*(const BezierSpline<parametric_dimension, PointTypeRHS,
                               ScalarRHS>& rhs)
      const {  // This multiplication operator is based on the algorithm
    // presented in the thesis from G. Elber (1992)

    // Initialize return value
    using PointTypeReturnT = decltype(PhysicalPointType{} * PointTypeRHS{});
    using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});

    // Determine the degrees of the resulting spline
    std::array<IndexingType, parametric_dimension> product_degrees;
    for (IndexingType param_dim{}; param_dim < parametric_dimension;
         param_dim++) {
      product_degrees[param_dim] = degrees[param_dim] + rhs.degrees[param_dim];
    }

    // Initialize the return type
    BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
        return_spline(product_degrees);

    // Start calculating the new control points
    for (IndexingType i{}; i < return_spline.NumberOfControlPoints; i++) {
      return_spline.product_combine_control_points(
          (*this), rhs, return_spline.local_to_global_index(i),
          static_cast<ScalarReturnT>(1.));
    }
    return return_spline;
  }
};

}  // namespace beziermanipulation

#endif  // SRC_BEZIER_HPP
