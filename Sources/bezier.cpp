#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

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

  template <typename... scalar>
  explicit constexpr Point(const scalar&... coords)
      : std::array<BaseType, space_dim>{coords...} {
    static_assert(sizeof...(coords) == space_dim,
                  "Base Logical Error: You are "
                  "instantiating a Point object with "
                  "more or less than "
                  "space_dim components.");
  }

  Point operator*(const double& scale) const {
    Point<space_dim, BaseType> product;
    for (unsigned int i = 0; i < space_dim; ++i)
      product.at(i) = (*this)[i] * scale;
    return product;
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

template <std::size_t parametric_dimension, std::size_t physical_dimension,
          typename Scalar>
class BezierSpline {
  using IndexingType = std::size_t;
  using PointTypePhysical = Point<physical_dimension, Scalar>;
  using PointTypeParametric = Point<parametric_dimension, Scalar>;

  constexpr void update_index_offsets() {
    index_offsets[0] = 1;
    for (unsigned int i{1}; i < parametric_dimension; i++)
      index_offsets[i] = index_offsets[i - 1] * (degrees[i - 1] + 1);
  }

 public:
  std::array<IndexingType, parametric_dimension> degrees;
  std::array<IndexingType, parametric_dimension> index_offsets;
  IndexingType n_ctps{};
  std::vector<PointTypePhysical> control_points{};

  constexpr BezierSpline( const BezierSpline &bezier_spline) = default;

  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg)
      : degrees{deg} {
    n_ctps = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++)
      n_ctps *= degrees[i] + 1;
    control_points.resize(n_ctps);
    update_index_offsets();
  };

  constexpr BezierSpline(
      const std::array<std::size_t, parametric_dimension> deg,
      const std::vector<PointTypePhysical> points)
      : degrees{deg}, control_points{points} {
    n_ctps = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++) n_ctps *= degrees[i] + 1;

    update_index_offsets();
    assert(n_ctps == points.size());
  };

  template <typename... T>
  constexpr PointTypePhysical& control_point(const T... index) {
    static_assert(sizeof...(T) == parametric_dimension,
                  "Unspecified number of indices.");
    unsigned int c_i{0}, i{};
    ((c_i += index_offsets[i++] * index), ...);
    return control_points[c_i];
  }

  constexpr BezierSpline& order_elevate_along_parametric_dimension(
      const IndexingType par_dim) {
    // Calculate index Offsets to facilitate working on 1D array
    const unsigned int n_starting_points = (n_ctps / (degrees[par_dim] + 1));
    const unsigned int starting_point_offset =
        index_offsets[par_dim] * (degrees[par_dim] + 1);
    const int starting_points_per_group = index_offsets[par_dim];
    const int n_groups = n_starting_points / starting_points_per_group;

    // Resize the CTPS vector accordingly
    n_ctps = n_ctps / (degrees[par_dim] + 1) * (degrees[par_dim] + 2);
    control_points.resize(n_ctps);
    degrees[par_dim]++;

    // Local Counter
    unsigned int global_index = n_ctps - 1;

    // Precalculations
    const Scalar inverse_factor =
        static_cast<Scalar>(1) / static_cast<Scalar>(degrees[par_dim]);
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
          const Scalar factor = static_cast<Scalar>(i) * inverse_factor;

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

  constexpr PointTypePhysical evaluate(
      const PointTypeParametric& par_coords) const {
    // Work on copy of control_point
    std::vector<PointTypePhysical> ctps_copy{control_points};
    IndexingType ctps_to_consider = n_ctps;

    for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
      Scalar factor = par_coords[par_dim];
      Scalar inv_factor = 1. - par_coords[par_dim];

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

  template <typename ScalarRHS>
  constexpr BezierSpline<parametric_dimension, physical_dimension,
                         decltype(Scalar{} * ScalarRHS{})>
  operator+(const BezierSpline<parametric_dimension, physical_dimension,
                               ScalarRHS>& rhs) {
    // Initialize return value
    using ScalarReturnT = decltype(Scalar{} * ScalarRHS{});
    BezierSpline<parametric_dimension, physical_dimension, ScalarReturnT>
        return_spline(degrees, control_points);

    bool rhs_needs_copy = false;
    for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
      rhs_needs_copy =
          rhs_needs_copy || (degrees[par_dim] > rhs.degrees[par_dim]);
    }
    if (rhs_needs_copy) {
      BezierSpline<parametric_dimension, physical_dimension, ScalarRHS>
          rhs_copy{rhs};
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
      for (IndexingType i_ctps{}; i_ctps < return_spline.n_ctps; i_ctps++){
        return_spline.control_points[i_ctps] += rhs.control_points[i_ctps];
      }
     return return_spline;
    }
  }
};

int main() {
  std::cout << "Working on the Bezier Splines\n";
  const unsigned int n_x{1}, n_y{1}, n_z{1};

  auto spline_test =
      BezierSpline<3, 3, double>(std::array<std::size_t, 3>{n_x, n_y, n_z});

  for (unsigned int i{0}; i <= n_x; i++) {
    for (unsigned int j{0}; j <= n_y; j++) {
      for (unsigned int k{0}; k <= n_z; k++) {
        spline_test.control_point(i, j, k) =
            Point<3, double>((double)i, (double)j, (double)k);
      }
    }
  }

  BezierSpline<3, 3, double> spline_test_copy {spline_test};

  spline_test.order_elevate_along_parametric_dimension(1);

  const auto sum_spline = spline_test + spline_test_copy;

  for (unsigned int i{0}; i < sum_spline.n_ctps; i++) {
    std::cout << "[";
    for (unsigned int dim{0}; dim < 3; dim++) {
      std::cout << " " << sum_spline.control_points[i][dim];
    }
    std::cout << " ] " << i << std::endl;
  }

  const auto eval_point = sum_spline.evaluate(Point<3, double>{0.2, .7, .5});
  for (unsigned int i{0}; i < 3; i++) std::cout << eval_point[i] << std::endl;

  return 0;
}