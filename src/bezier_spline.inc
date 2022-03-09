template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr void BezierSpline<parametric_dimension, PhysicalPointType,
                            ScalarType>::update_index_offsets() {
  index_offsets[0] = 1;
  for (unsigned int i{1}; i < parametric_dimension; i++)
    index_offsets[i] = index_offsets[i - 1] * (degrees[i - 1] + 1);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::size_t, parametric_dimension> BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::local_to_global_index(const IndexingType& local_index) const {
  std::array<IndexingType, parametric_dimension> indexList{};
  for (unsigned int i{0}; i < parametric_dimension; i++) {
    indexList[i] = (local_index / index_offsets[i]) % (degrees[i] + 1);
  }
  return indexList;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PhysicalPointLHS, typename ScalarLHS,
          typename PhysicalPointRHS, typename ScalarRHS, typename... T>
constexpr void
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    product_combine_control_points(
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
          static_cast<IndexingType>(indices)..., static_cast<IndexingType>(i)};
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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::control_point(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  unsigned int c_i{0}, i{};
  ((c_i += index_offsets[i++] * index), ...);
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::control_point(const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  unsigned int c_i{0}, i{};
  ((c_i += index_offsets[i++] * index), ...);
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    control_point(
        const std::array<IndexingType, parametric_dimension>& index) const {
  unsigned int c_i{0}, i{};
  for (unsigned int i{}; i < parametric_dimension; i++) {
    c_i += index_offsets[i] * index[i];
  }
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    control_point(const std::array<IndexingType, parametric_dimension>& index) {
  unsigned int c_i{0}, i{};
  for (unsigned int i{}; i < parametric_dimension; i++) {
    c_i += index_offsets[i] * index[i];
  }
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    order_elevate_along_parametric_dimension(const IndexingType par_dim) {
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
      static_cast<ScalarType_>(1) / static_cast<ScalarType_>(degrees[par_dim]);
  const IndexingType variable_offset_factor =
      index_offsets[par_dim] * (degrees[par_dim] - 1);

  // Vector is calculated from back to front, to hinder overwrite
  for (int group_index{n_groups - 1}; group_index >= 0; group_index--) {
    // Local variables
    const unsigned int first_index_in_group =
        group_index * (starting_point_offset);
    IndexingType i = degrees[par_dim] - 1;

    // Fix the last entry for element in the group
    for (int index_in_group{starting_points_per_group - 1}; index_in_group >= 0;
         index_in_group--) {
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
        const ScalarType_ factor = static_cast<ScalarType_>(i) * inverse_factor;

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
    for (int index_in_group{starting_points_per_group - 1}; index_in_group >= 0;
         index_in_group--) {
      // No we actually start the algorithm
      control_points[global_index] =
          control_points[first_index_in_group + index_in_group];
      global_index--;
    }
  }
  update_index_offsets();
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::evaluate(const PointTypeParametric_& par_coords)
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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator==(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        rhs) const {
  // Check if degrees fit and if control_points are the same
  if (rhs.degrees == degrees) {
    return rhs.control_points == control_points;
  } else {
    return false;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*=(
    const ScalarType& scalar) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] *= scalar;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+=(
    const PhysicalPointType& point_shift) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] += point_shift;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+(
    const PhysicalPointType& point_shift) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline += point_shift;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-=(
    const PhysicalPointType& point_shift) {
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    control_points[i] -= point_shift;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-(
    const PhysicalPointType& point_shift) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline -= point_shift;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*(
    const ScalarType& scalar) const {
  BezierSpline scaled_spline{(*this)};
  scaled_spline *= scalar;
  return scaled_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr BezierSpline<parametric_dimension,
                       decltype(PhysicalPointType{} + PointTypeRHS{}),
                       decltype(ScalarType{} + ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+(
    const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs)
    const {
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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::ExtractDimension(unsigned int dimension) const {
  assert(dimension < PointTypePhysical_::kSpatialDimension);
  BezierSpline<parametric_dimension, ScalarType, ScalarType> extracted_spline(
      degrees);
  for (std::size_t i{}; i < NumberOfControlPoints; i++) {
    extracted_spline.control_points[i] = control_points[i][dimension];
  }
  return extracted_spline;
}
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::power(
    const unsigned int power) const {
  // @TODO use log2(power) algorithm that squares the result to minimize
  // multiplications
  static_assert(std::is_scalar_v<PhysicalPointType>,
                "Only Scalar-type Splines can be raised to a power.");
  assert(("Not implemented, as raising to 0 would be inefficient.", power > 0));
  BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
      power_spline{(*this)};
  for (int i{1}; i < power; i++) {
    power_spline = power_spline * (*this);
  }
  return power_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr BezierSpline<parametric_dimension,
                       decltype(PhysicalPointType{} * PointTypeRHS{}),
                       decltype(ScalarType{} * ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*(
    const BezierSpline<parametric_dimension, PointTypeRHS,
                       ScalarRHS>& rhs)
    const {  // This multiplication operator is based on the algorithm
  // presented in the thesis from G. Elber (1992)

  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} * PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType{} * ScalarRHS{});

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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr BezierSpline<parametric_dimension_inner_spline, PhysicalPointType,
                       decltype(ScalarType{} * ScalarRHS{})>

BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::compose(
    const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
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

  // Extract splines only once, and store data on heap. The size corresponds
  // to the parametric dimensions of the outer spline representation
  std::vector<BezierSpline<parametric_dimension_inner_spline, ScalarReturnT,
                           ScalarReturnT>>
      inner_spline_coordinate_seperations{parametric_dimension};
  for (std::size_t i_outer_parametric_dimension{};
       i_outer_parametric_dimension < parametric_dimension;
       i_outer_parametric_dimension++) {
    // Extract the current dimension of the spline
    inner_spline_coordinate_seperations[i_outer_parametric_dimension] =
        inner_function.ExtractDimension(i_outer_parametric_dimension);
  }  // Val_Jac

  // Loop over control points of the outer spline to start the interpolation
  for (unsigned int i_control_point_outer_spline{};
       i_control_point_outer_spline < NumberOfControlPoints;
       i_control_point_outer_spline++) {
    // Retrieve indices
    const auto outer_spline_ctps_indices =
        local_to_global_index(i_control_point_outer_spline);

    // Create a 2D vector of all control-point factors to be filled later on
    std::vector<BezierSpline<parametric_dimension_inner_spline, ScalarReturnT,
                             ScalarReturnT>>
        factorization_splines{parametric_dimension};

    // Fill the factor list
    for (std::size_t i_outer_parametric_dimension{};
         i_outer_parametric_dimension < parametric_dimension;
         i_outer_parametric_dimension++) {
      const auto i_basis_function_index =
          outer_spline_ctps_indices[i_outer_parametric_dimension];
      const auto n_outer_degree =
          degrees[i_outer_parametric_dimension];  // Reference for readability
      const auto& inner_spline_xi =
          inner_spline_coordinate_seperations[i_outer_parametric_dimension];

      // Boundary cases
      if (i_basis_function_index == 0) {
        factorization_splines[i_outer_parametric_dimension] =
            ((1 - inner_spline_xi).power(n_outer_degree));
      } else if (i_basis_function_index == n_outer_degree) {
        factorization_splines[i_outer_parametric_dimension] =
            (inner_spline_xi.power(n_outer_degree));
      } else {
        factorization_splines[i_outer_parametric_dimension] =
            (utils::FastBinomialCoefficient::choose(n_outer_degree,
                                                    i_basis_function_index) *
             inner_spline_xi.power(i_basis_function_index) *
             (1 - inner_spline_xi)
                 .power(n_outer_degree - i_basis_function_index));
      }
    }

    // Multiply all factor splines with each other
    BezierSpline<parametric_dimension_inner_spline, ScalarReturnT,
                 ScalarReturnT>
        control_point_multiplication_factor{factorization_splines[0]};
    for (std::size_t i_outer_parametric_dimension{1};
         i_outer_parametric_dimension < parametric_dimension;
         i_outer_parametric_dimension++) {
      control_point_multiplication_factor =
          control_point_multiplication_factor *
          factorization_splines[i_outer_parametric_dimension];
    }

    // Now that all factors have been precomputed we can determine the
    // control-points of the resulting spline
    for (IndexingType i_control_point_composition{};
         i_control_point_composition < composition.NumberOfControlPoints;
         i_control_point_composition++) {
      // Start interpolation
      composition.control_points[i_control_point_composition] +=
          control_point_multiplication_factor
              .control_points[i_control_point_composition] *
          control_points[i_control_point_outer_spline];
    }
  }

  return composition;
}