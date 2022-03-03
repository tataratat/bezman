#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

// Testing
#include <random>



template <typename Scalar, unsigned int max_degree>
class BinomialCoefficientLookupTableCreator {
 public:
  static constexpr std::size_t n_values = max_degree * (max_degree + 1) / 2;

  static constexpr std::array<Scalar, n_values> calculate_table() {
    std::array<Scalar, n_values> ret_val{};
    ret_val[0] = static_cast<Scalar>(1);
    for (std::size_t i{0}; i < max_degree; i++) {
      std::size_t offset = i * (i + 1) / 2;
      for (std::size_t j{0}; j <= i; j++) {
        if (i == j || j == 0) {
          ret_val[offset + j] = static_cast<Scalar>(1);
        } else {
          ret_val[offset + j] =
              ret_val[offset - i + j - 1] + ret_val[offset - i + j];
        }
      }
    }
    return ret_val;
  }
};

// MINIMALISTIC COPY FROM CAMPIGA POINTS
template <unsigned int space_dim, typename BaseType = double>
class Point : public std::array<BaseType, space_dim> {
 public:
  constexpr Point(const Point&) = default;

  constexpr Point() = default;

  template <typename... scalar>
  explicit constexpr Point(const scalar&... coords)
      : std::array<BaseType, space_dim>{coords...} {
    static_assert(sizeof...(coords) == space_dim,
                  "Base Logical Error: You are "
                  "instantiating a Point object with "
                  "more or less than "
                  "space_dim components.");
  }

  constexpr Point operator+(const Point& rhs) const {
    Point<space_dim, BaseType> sum{(*this)};
    sum += rhs;
    return sum;
  }

  constexpr Point& operator+=(const Point& rhs) {
    for (unsigned int i = 0; i < space_dim; i++) {
      (*this)[i] = (*this)[i] + rhs[i];
    }
    return (*this);
  }

  constexpr void scale(const double& scale) {
    for (unsigned int i = 0; i < space_dim; ++i) (*this)[i] *= scale;
  }

  Point operator*(const double& scale) const {
    Point<space_dim, BaseType> product;
    for (unsigned int i = 0; i < space_dim; ++i)
      product.at(i) = (*this)[i] * scale;
    return product;
  }

  template<typename Scalar>
  decltype(Scalar{} * BaseType{}) 
  operator*(const Point<space_dim,Scalar>& point) const {
    using ScalarReturn = decltype(Scalar{} * BaseType{});
    ScalarReturn result{};
    for (unsigned int i{}; i < space_dim; i++){
      result += (*this)[i] * point[i];
    }
    return result;
  }

  friend std::ostream& operator<<(std::ostream& os, const Point& p) {
    for (size_t i = 0; i < space_dim; ++i) {
      os << (i == 0 ? "(" : ", ");
      os << std::setw(5) << std::setprecision(3) << p[i];
      os << (i == space_dim - 1 ? ")" : "");
    }
    return os;
  }
};

#ifndef MAX_BINOMIAL_DEGREE
#define MAX_BINOMIAL_DEGREE 30u
#endif

// Global Lookup
class FastBinomialCoefficient {
 private:
  using ScalarType = std::size_t;
  constexpr static auto n_values =
      BinomialCoefficientLookupTableCreator<ScalarType,
                                            MAX_BINOMIAL_DEGREE>::n_values;
  constexpr static std::array<ScalarType, n_values> look_up =
      BinomialCoefficientLookupTableCreator<
          ScalarType, MAX_BINOMIAL_DEGREE>::calculate_table();

 public:
  constexpr static ScalarType choose(const unsigned int n,
                                     const unsigned int i) {
    assert(n < MAX_BINOMIAL_DEGREE);
    return look_up[n * (n + 1) / 2 + i];
  }
};

template< std::size_t parametric_dimension, 
          typename PhysicalPointType,
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
   local_to_global_index(const IndexingType &local_index)
   {
     std::array<IndexingType, parametric_dimension> indexList{};
     for (unsigned int i{0}; i < parametric_dimension; i++)
     {
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
    for (int i{std::max(0, k - n)};
         i <= std::min(k, m); i++) {
      // Calculate Factor
      const ScalarType lFactor =
          static_cast<ScalarType>(FastBinomialCoefficient::choose(m, i) *
                                  FastBinomialCoefficient::choose(n, k - i)) /
          static_cast<ScalarType>(FastBinomialCoefficient::choose(m + n, k));

      // Now decide if continue recursion
      if constexpr ((depth + 1) == parametric_dimension){
          const std::array<IndexingType, parametric_dimension> ind_lhs{static_cast<IndexingType>(indices)..., static_cast<IndexingType>(i)};
          std::array<IndexingType,
                     parametric_dimension> ind_rhs{};
          for (unsigned int j{}; j < parametric_dimension; j++) {
            ind_rhs[j] = ctpsIndex[j] - ind_lhs[j];
          }
          (*this).control_point(ctpsIndex) += P_spline.control_point(ind_lhs) *
                                          Q_spline.control_point(ind_rhs) *
                                          factor * lFactor;
      } else {
        product_combine_control_points(P_spline, Q_spline, ctpsIndex, factor * lFactor, indices..., i);
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
    const unsigned int n_starting_points = (NumberOfControlPoints / (degrees[par_dim] + 1));
    const unsigned int starting_point_offset =
        index_offsets[par_dim] * (degrees[par_dim] + 1);
    const int starting_points_per_group = index_offsets[par_dim];
    const int n_groups = n_starting_points / starting_points_per_group;

    // Resize the CTPS vector accordingly
    NumberOfControlPoints = NumberOfControlPoints / (degrees[par_dim] + 1) * (degrees[par_dim] + 2);
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

  /// Evaluate the spline via the deCasteljau algorithm
  constexpr PointTypePhysical_ evaluate(
      const PointTypeParametric_& par_coords) const {
    // Work on copy of control_point
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
  operator+(
      const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs) {
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

    // If the RHS needs a copy call function recursively
    if (rhs_needs_copy) {
      BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS> rhs_copy{rhs};
      for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
        while (degrees[par_dim] > rhs_copy.degrees[par_dim]) {
          rhs_copy.order_elevate_along_parametric_dimension(par_dim);
        }
      }
      return (*this) + rhs_copy;
    } else {
      for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
        while (rhs.degrees[par_dim] < return_spline.degrees[par_dim]) {
          return_spline.order_elevate_along_parametric_dimension(par_dim);
        }
      }
      for (IndexingType i_ctps{}; i_ctps < return_spline.NumberOfControlPoints; i_ctps++) {
        return_spline.control_points[i_ctps] += rhs.control_points[i_ctps];
      }
      return return_spline;
    }
  }

  template <typename PointTypeRHS, typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension,
                         decltype(PhysicalPointType{} * PointTypeRHS{}),
                         decltype(ScalarType_{} * ScalarRHS{})>
  operator*(
      const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs) {
    // This multiplication operator is based on the algorithm
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
    for (IndexingType i{}; i < return_spline.NumberOfControlPoints; i++){
      return_spline.product_combine_control_points(
        (*this), 
        rhs, 
        return_spline.local_to_global_index(i), 
        static_cast<ScalarReturnT>(1.));
    }
    return return_spline;
  }
};

int main() {
  std::cout << "Working on the Bezier Splines\n";
  const unsigned int n_x{1}, n_y{1}, n_z{1};

  auto spline_test =
      BezierSpline<3u, Point<3,double>, double>(std::array<std::size_t, 3>{n_x, n_y, n_z});

  for (unsigned int i{0}; i <= n_x; i++) {
    for (unsigned int j{0}; j <= n_y; j++) {
      for (unsigned int k{0}; k <= n_z; k++) {
        spline_test.control_point(i, j, k) =
            Point<3, double>((double)i, (double)j, (double)k);
      }
    }
  }

  BezierSpline<3u, Point<3, double>, double> spline_test_copy{spline_test};

  const auto sum_spline = spline_test * spline_test_copy;

  for (unsigned int i{0}; i < sum_spline.NumberOfControlPoints; i++) {
      std::cout << sum_spline.control_points[i];
    std::cout << " : " << i << std::endl;
  }

  for (int i{}; i < 10; i++){
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        z{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto eval_point = sum_spline.evaluate(Point<3, double>{x,y,z});
    const auto eval_point1 = spline_test.evaluate(Point<3, double>{x, y, z});
    const auto eval_point2 =
        spline_test_copy.evaluate(Point<3, double>{x, y, z});

    std::cout << std::setw(5) << std::setprecision(3) 
              << "Result :\t"
              << eval_point << "\tVec1 :\t" << eval_point1 << "\tVec2 :\t"
              << eval_point2 << "\tExpected Result :\t"
              << eval_point1 * eval_point2 << std::endl;
  }

  // // Second Test
  // BezierSpline<1u, Point<3, double>, double> line0{{1}};
  // BezierSpline<1u, Point<3, double>, double> line1{{2}};

  // line0.control_point(0) = Point<3u, double>{1, 0, 0};
  // line0.control_point(1) = Point<3u, double>{1, 1, 0};
  // line1.control_point(0) = Point<3u, double>{1, 0, 0};
  // line1.control_point(1) = Point<3u, double>{1, 1, 0};
  // line1.control_point(2) = Point<3u, double>{1, 2, 0};

  // const auto line_product = line0 * line1;

  // for (unsigned int i{0}; i < line_product.NumberOfControlPoints; i++) {
  //   std::cout << line_product.control_points[i];
  //   std::cout << " : " << i << std::endl;
  // }

  // for (int i{}; i < 10; i++) {
  //   const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  //   const auto eval_point = line_product.evaluate(Point<1, double>{x});
  //   const auto eval_point1 = line0.evaluate(Point<1, double>{x});
  //   const auto eval_point2 = line1.evaluate(Point<1, double>{x});

  //   std::cout << std::setw(5) << std::setprecision(3) << "Result :\t"
  //             << eval_point << "\tVec1 :\t" << eval_point1 << "\tVec2 :\t"
  //             << eval_point2 << "\tExpected Result :\t"
  //             << eval_point1 * eval_point2 << std::endl;
  // }

  return 0;
}