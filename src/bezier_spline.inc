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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr void BezierSpline<parametric_dimension, PhysicalPointType,
                            ScalarType>::UpdateIndexOffsets_() {
  index_offsets[0] = 1;
  for (unsigned int i{1}; i < parametric_dimension; i++)
    index_offsets[i] = index_offsets[i - 1] * (degrees[i - 1] + 1);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::size_t, parametric_dimension> BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::LocalToGlobalIndex(const IndexingType& local_index) const {
  std::array<IndexingType, parametric_dimension> indexList{};
  for (unsigned int i{0}; i < parametric_dimension; i++) {
    indexList[i] = (local_index / index_offsets[i]) % (degrees[i] + 1);
  }
  return indexList;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr typename BezierSpline<parametric_dimension, PhysicalPointType,
                                ScalarType>::IndexingType_
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    GlobalToLocalIndex(const std::array<std::size_t, parametric_dimension>&
                           local_index) const {
  IndexingType c_i{local_index[0]};
  for (unsigned int i{1}; i < parametric_dimension; i++) {
    c_i += index_offsets[i] * local_index[i];
  }
  return c_i;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PhysicalPointLHS, typename ScalarLHS,
          typename PhysicalPointRHS, typename ScalarRHS, typename... T>
constexpr void
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    CombineControlPointsForProduct_(
        const BezierSpline<parametric_dimension, PhysicalPointLHS, ScalarLHS>&
            P_spline,
        const BezierSpline<parametric_dimension, PhysicalPointRHS, ScalarRHS>&
            Q_spline,
        const std::array<IndexingType, parametric_dimension>& ctpsIndex,
        const ScalarType factor, const T&... indices) {
  // Some constant indices and degrees
  const int depth = sizeof...(indices);
  const int k = ctpsIndex[depth];
  const int m = P_spline.GetDegrees()[depth];
  const int n = Q_spline.GetDegrees()[depth];

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
      (*this).ControlPoint(ctpsIndex) += P_spline.ControlPoint(ind_lhs) *
                                         Q_spline.ControlPoint(ind_rhs) *
                                         factor * lFactor;
    } else {
      CombineControlPointsForProduct_(P_spline, Q_spline, ctpsIndex,
                                      factor * lFactor, indices..., i);
    }
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const T... index) const {
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
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  unsigned int c_i{0}, i{};
  ((c_i += index_offsets[i++] * index), ...);
  return control_points[c_i];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const std::array<IndexingType, parametric_dimension>& index) const {
  return control_points[GlobalToLocalIndex(index)];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::ControlPoint(
    const std::array<IndexingType, parametric_dimension>& index) {
  return control_points[GlobalToLocalIndex(index)];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    OrderElevateAlongParametricDimension(const IndexingType par_dim) {
  // Calculate index Offsets to facilitate working on 1D array
  const unsigned int n_starting_points =
      (GetNumberOfControlPoints() / (degrees[par_dim] + 1));
  const unsigned int starting_point_offset =
      index_offsets[par_dim] * (degrees[par_dim] + 1);
  const int starting_points_per_group = index_offsets[par_dim];
  const int n_groups = n_starting_points / starting_points_per_group;

  // Resize the CTPS vector accordingly
  number_of_control_points = number_of_control_points / (degrees[par_dim] + 1) *
                             (degrees[par_dim] + 2);
  control_points.resize(GetNumberOfControlPoints());
  degrees[par_dim]++;

  // Local Counter
  unsigned int global_index = GetNumberOfControlPoints() - 1;

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
  UpdateIndexOffsets_();
  utils::Logger::ExtendedInformation(
      "Spline elevating along parametric dimension " + std::to_string(par_dim) +
      " to degree " + std::to_string(degrees[par_dim]));
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    DerivativeWRTParametricDimension(
        const std::array<IndexingType, parametric_dimension>& orders) const {
  utils::Logger::Logging(
      "Calculating derivative splines with orders " +
          [](const std::array<IndexingType, parametric_dimension>& orders)
          -> std::string {
        std::string result("[ " + std::to_string(orders[0]));
        for (IndexingType i{1}; i < parametric_dimension; i++) {
          result += ", " + std::to_string(orders[i]);
        }
        return result + "]";
      }(orders));

  // Create a copy of the new degrees which will be updated in every iteration
  auto new_degrees = degrees;
  // Create a copy of control points (later be used for in place operations)
  auto new_control_points = GetControlPoints();

  // Determine how many control points are expected in the current spline
  auto compute_number_of_control_points = [&new_degrees]() -> IndexingType {
    IndexingType number_of_control_points_ = 1u;
    for (unsigned int i{}; i < parametric_dimension; i++) {
      number_of_control_points_ *= new_degrees[i] + 1;
    }
    return number_of_control_points_;
  };

  // Determine control point IDs using current values of new_degrees
  auto global_id_with_new_degrees = [&new_degrees](IndexingType global_id)
      -> std::array<IndexingType, parametric_dimension> {
    std::array<IndexingType, parametric_dimension> multi_index{
        global_id % (new_degrees[0] + 1)};
    for (unsigned int i_para{1}; i_para < parametric_dimension; i_para++) {
      global_id /= (new_degrees[i_para - 1] + 1);
      multi_index[i_para] = global_id % (new_degrees[i_para] + 1);
    }
    return multi_index;
  };

  // Loop over order vector and apply
  for (IndexingType i_para_dim{}; i_para_dim < parametric_dimension;
       i_para_dim++) {
    if (orders[i_para_dim] == 0) {
      continue;
    }
    // Update the new degrees
    new_degrees[i_para_dim] -= orders[i_para_dim];

    // Compute the constant factor p!/k!
    const ScalarType& factor = [](const IndexingType& p,
                                  const IndexingType& k) -> ScalarType {
      ScalarType factor{static_cast<ScalarType>(p)};
      for (IndexingType i{1}; i < k; i++) {
        factor *= p - i;
      }
      return factor;
    }(degrees[i_para_dim], orders[i_para_dim]);

    // Start looping over all control points in derivative and determine
    // their new values
    for (IndexingType i_point{0}; i_point < compute_number_of_control_points();
         i_point++) {
      auto multi_index = global_id_with_new_degrees(i_point);
      const auto i_point_index = GlobalToLocalIndex(multi_index);
      for (IndexingType j_point_offset{1}; j_point_offset <= orders[i_para_dim];
           j_point_offset++) {
        multi_index[i_para_dim]++;
        const ScalarType local_coefficient =
            (j_point_offset % 2 == 0 ? 1.0 : -1.0) *
            utils::FastBinomialCoefficient::choose(orders[i_para_dim],
                                                   j_point_offset);
        new_control_points[i_point_index] +=
            local_coefficient *
            new_control_points[GlobalToLocalIndex(multi_index)];
      }
      new_control_points[i_point_index] *=
          (orders[i_para_dim] % 2 == 0 ? 1.0 : -1.0) * factor;
    }
  }

  // Create new spline and
  BezierSpline derivative{new_degrees};

  // Copy all control points into the vector
  for (IndexingType i_point{}; i_point < derivative.GetNumberOfControlPoints();
       i_point++) {
    derivative.GetControlPoints()[i_point] =
        new_control_points[GlobalToLocalIndex(
            derivative.LocalToGlobalIndex(i_point))];
  }

  utils::Logger::Logging("Derivative successfully computed");
  return derivative;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Evaluate(
    const PointTypeParametric_& par_coords) const {
  /*
   * Notation:
   * A group the number of all starting points that are subsequent, i.e. the
   * size of a group is the number of subsequent 0s. E.g. in the example:
   * Group size of index 0 is 1, group size of index 1 is 3, etc.. Group
   * offset is thus the distance between two group-starting-points. The idea is
   * to collapse the groups one by one starting from the first parametric
   * dimension going forward. For indexing, see SplitAtPosition(...)
   */

  // Work on copy of ControlPoint
  std::vector<PhysicalPointType_> ctps_copy{control_points};
  IndexingType ctps_to_consider = GetNumberOfControlPoints();

  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    ScalarType_ factor = par_coords[par_dim];
    ScalarType_ inv_factor = 1. - par_coords[par_dim];

    ctps_to_consider /= degrees[par_dim] + 1;

    // For every starting position
    for (IndexingType starting_point_id{0};
         starting_point_id < ctps_to_consider; starting_point_id++) {
      const auto offset =
          (index_offsets[par_dim] * (degrees[par_dim] + 1)) * starting_point_id;
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
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    BasisFunctionContributions(const PointTypeParametric_& par_coords) const {
  std::array<std::vector<ScalarType>, parametric_dimension> factors;

  // Preallocate space to vectors
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    factors[par_dim] =
        utils::algorithms::BernsteinPolynomial<ScalarType>::Evaluate(
            degrees[par_dim], par_coords[par_dim]);
  }

  return factors;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::vector<ScalarType> BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::BasisFunctions(const PointTypeParametric_& par_coords) const {
  // Recursively combine basis function values per parametric dimension
  return utils::algorithms::RecursiveCombine(
      BasisFunctionContributions(par_coords));
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType BezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::ForwardEvaluate(const PointTypeParametric_& par_coords) const {
  /*
   * Here we avoid Point operations and try to precompute as many values as
   * possible, this greatly decreases the complexity of the necessary code.
   */
  const auto factors = BasisFunctions(par_coords);

  /*
   * Now that all values have been precomputed we have to sum up the
   * individual contributions using (compile-time) recursion
   */
  PhysicalPointType_ evaluation_point{};

  for (IndexingType i_basis_function{0}; i_basis_function < factors.size();
       i_basis_function++) {
    evaluation_point +=
        factors[i_basis_function] * GetControlPoints()[i_basis_function];
  }

  return evaluation_point;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    BasisFunctionsDerivativeContributions(
        const PointTypeParametric_& par_coords,
        const std::array<std::size_t, parametric_dimension>& nth_derivs) const {
  std::array<std::vector<ScalarType>, parametric_dimension> factors;

  // Preallocate space to vectors
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    factors[par_dim] =
        utils::algorithms::BernsteinPolynomial<ScalarType>::EvaluateDerivative(
            degrees[par_dim], par_coords[par_dim], nth_derivs[par_dim]);
  }

  return factors;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::vector<ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    BasisFunctionsDerivatives(
        const PointTypeParametric_& par_coords,
        const std::array<std::size_t, parametric_dimension>& nth_derivs) const {
  std::array<std::vector<ScalarType>, parametric_dimension> factors;

  // Preallocate space to vectors
  return utils::algorithms::RecursiveCombine(
      BasisFunctionsDerivativeContributions(par_coords, nth_derivs));
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    EvaluateDerivative(
        const PointTypeParametric_& par_coords,
        const std::array<std::size_t, parametric_dimension>& nth_derivs) const {
  /*
   * Here we avoid Point operations and try to precompute as many values as
   * possible, this greatly decreases the complexity of the necessary code.
   */
  const auto factors = BasisFunctionsDerivatives(par_coords, nth_derivs);

  /*
   * Now that all values have been precomputed we have to sum up the
   * individual contributions using (compile-time) recursion
   */
  PhysicalPointType_ evaluation_point{};

  for (IndexingType i_basis_function{0}; i_basis_function < factors.size();
       i_basis_function++) {
    evaluation_point +=
        factors[i_basis_function] * GetControlPoints()[i_basis_function];
  }

  return evaluation_point;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator==(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        rhs) const {
  // Check if degrees fit and if control_points are the same
  if (rhs.GetDegrees() == degrees) {
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
  for (std::size_t i{}; i < GetNumberOfControlPoints(); i++) {
    control_points[i] *= scalar;
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+=(
    const PhysicalPointType& point_shift) {
  for (std::size_t i{}; i < GetNumberOfControlPoints(); i++) {
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
  for (std::size_t i{}; i < GetNumberOfControlPoints(); i++) {
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
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-()
    const {
  BezierSpline inverted_spline{(*this)};
  for (IndexingType i{}; i < GetNumberOfControlPoints(); i++) {
    inverted_spline.control_points[i] = -control_points[i];
  }
  return inverted_spline;
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
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+(
    const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs)
    const {
  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} + PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});

  // Return type
  BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
      return_spline(degrees, control_points);

  // Check if the right hand side requires a copy as it should not be
  // altered for this purpose
  bool rhs_needs_copy = false;
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    rhs_needs_copy =
        rhs_needs_copy || (degrees[par_dim] > rhs.GetDegrees()[par_dim]);
  }

  // Increase the order of the copied spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > return_spline.GetDegrees()[par_dim]) {
      return_spline.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  if (!rhs_needs_copy) {
    for (IndexingType i_ctps{};
         i_ctps < return_spline.GetNumberOfControlPoints(); i_ctps++) {
      return_spline.control_points[i_ctps] += rhs.control_points[i_ctps];
    }
    return return_spline;
  } else {
    // use commutativity of addition to create a copy of rhs
    return rhs + return_spline;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator+=(
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType> rhs) {
  // Increase the order of (*this) spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > degrees[par_dim]) {
      OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Do the same thing to the rhs
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] < degrees[par_dim]) {
      rhs.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Add control points
  for (IndexingType i_ctps{}; i_ctps < GetNumberOfControlPoints(); i_ctps++) {
    control_points[i_ctps] += rhs.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-(
    BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS> rhs) const {
  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} + PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType_{} * ScalarRHS{});

  // Return type
  BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
      return_spline(degrees, control_points);

  // Increase the order of return_spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > degrees[par_dim]) {
      return_spline.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Do the same thing to the rhs
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] < degrees[par_dim]) {
      rhs.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Add control points
  for (IndexingType i_ctps{}; i_ctps < GetNumberOfControlPoints(); i_ctps++) {
    return_spline.control_points[i_ctps] -= rhs.control_points[i_ctps];
  }
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator-=(
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType> rhs) {
  // Increase the order of (*this) spline to be greater or equal to the
  // RHSs order
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] > degrees[par_dim]) {
      OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Do the same thing to the rhs
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    while (rhs.GetDegrees()[par_dim] < degrees[par_dim]) {
      rhs.OrderElevateAlongParametricDimension(par_dim);
    }
  }
  // Add control points
  for (IndexingType i_ctps{}; i_ctps < GetNumberOfControlPoints(); i_ctps++) {
    control_points[i_ctps] -= rhs.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::ExtractDimension(const IndexingType& dimension)
    const {
  if constexpr (std::is_arithmetic_v<PhysicalPointType>) {
    if (dimension != 0ul) {
      utils::Logger::TerminatingError(
          "Extraction of Scalar-value spline requested");
    }
    return (*this);
  } else {
    if (dimension >= PhysicalPointType_::kSpatialDimension) {
      utils::Logger::TerminatingError(
          "Dimension out of range for Scalar extraction");
    }
    BezierSpline<parametric_dimension, ScalarType, ScalarType> extracted_spline(
        degrees);
    for (std::size_t i{}; i < GetNumberOfControlPoints(); i++) {
      extracted_spline.control_points[i] = control_points[i][dimension];
    }
    return extracted_spline;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::RaisePower(
    const IndexingType power) const {
  // Implements Recursive function of integer powers (can be storage intense)
  if constexpr (std::is_scalar_v<PhysicalPointType>) {
    if (power == 0) {
      utils::Logger::Error("Raising to 0 is inefficient. Please revisit");
      return BezierSpline<parametric_dimension, ScalarType_, ScalarType_>{
          // Degree 0 (use zero default initialization)
          std::array<std::size_t, parametric_dimension>{},
          // Set the single remaining ctps to 1
          std::vector<ScalarType>{1.}};
    }
    if (power == 1) {
      return BezierSpline<parametric_dimension, ScalarType_, ScalarType_>{
          (*this)};
    }

    // Start Recursion
    const auto recursive_power =
        RaisePower(static_cast<IndexingType>(power / 2));
    if ((power % 2) == 0) {
      return recursive_power * recursive_power;
    } else {
      return recursive_power * recursive_power * (*this);
    }
  } else {
    utils::Logger::TerminatingError(
        "Only Scalar-type Splines can be raised to a power.");
    return BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>{};
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::MaximumCorner() const {
  if (GetNumberOfControlPoints() == 0) {
    utils::Logger::TerminatingError(
        "Spline is non-defined or has no control-points");
  }
  PhysicalPointType maximum{control_points[0]};
  for (IndexingType i{1}; i < GetNumberOfControlPoints(); i++) {
    if constexpr (std::is_arithmetic_v<PhysicalPointType>) {
      maximum = std::max(maximum, control_points[i]);
    } else {
      for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
           i_dim++) {
        maximum[i_dim] = std::max(maximum[i_dim], control_points[i][i_dim]);
      }
    }
  }
  return maximum;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::MinimumCorner() const {
  if (GetNumberOfControlPoints() == 0) {
    utils::Logger::TerminatingError(
        "Spline is non-defined or has no control-points");
  }
  PhysicalPointType minimum{control_points[0]};
  for (IndexingType i{1}; i < GetNumberOfControlPoints(); i++) {
    if constexpr (std::is_arithmetic_v<PhysicalPointType>) {
      minimum = std::min(minimum, control_points[i]);
    } else {
      for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
           i_dim++) {
        minimum[i_dim] = std::min(minimum[i_dim], control_points[i][i_dim]);
      }
    }
  }
  return minimum;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    TransposeAndScale(const PhysicalPointType& transposition,
                      const PhysicalPointType& scale_vector) {
  utils::Logger::ExtendedInformation(
      "Starting scaling and transposing spline to move bounds");
  for (IndexingType i{}; i < GetNumberOfControlPoints(); i++) {
    if constexpr (std::is_scalar_v<PhysicalPointType>) {
      control_points[i] = (transposition + control_points[i]) * scale_vector;
    } else {
      for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
           i_dim++) {
        control_points[i][i_dim] =
            (transposition[i_dim] + control_points[i][i_dim]) *
            scale_vector[i_dim];
      }
    }
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::FitIntoUnitCube() {
  const auto maximum_corner = MaximumCorner();
  const auto minimum_corner = MinimumCorner();
  auto scale_vector = maximum_corner - minimum_corner;
  if constexpr (std::is_scalar_v<PhysicalPointType>) {
    scale_vector = static_cast<ScalarType>(1) / scale_vector;
  } else {
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      scale_vector[i_dim] = static_cast<ScalarType>(1) / scale_vector[i_dim];
    }
  }
  return TransposeAndScale(-minimum_corner, scale_vector);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool BezierSpline<parametric_dimension, PhysicalPointType,
                            ScalarType>::FitsIntoUnitCube() const {
  const auto maximum_corner = MaximumCorner();
  const auto minimum_corner = MinimumCorner();
  if constexpr (std::is_arithmetic_v<PhysicalPointType>) {
    if ((maximum_corner > static_cast<ScalarType>(1)) ||
        (minimum_corner < static_cast<ScalarType>(0))) {
      return false;
    }
  } else {
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      if ((maximum_corner[i_dim] > static_cast<ScalarType>(1)) ||
          (minimum_corner[i_dim] < static_cast<ScalarType>(0))) {
        return false;
      }
    }
  }
  return true;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::operator*(
    const BezierSpline<parametric_dimension, PointTypeRHS, ScalarRHS>& rhs)
    const {
  // This multiplication operator is based on the algorithm
  // presented in the thesis from G. Elber (1992)

  // Initialize return value
  using PointTypeReturnT = decltype(PhysicalPointType{} * PointTypeRHS{});
  using ScalarReturnT = decltype(ScalarType{} * ScalarRHS{});

  // Determine the degrees of the resulting spline
  std::array<IndexingType, parametric_dimension> product_degrees;
  for (IndexingType param_dim{}; param_dim < parametric_dimension;
       param_dim++) {
    product_degrees[param_dim] =
        degrees[param_dim] + rhs.GetDegrees()[param_dim];
  }

  // Initialize the return type
  BezierSpline<parametric_dimension, PointTypeReturnT, ScalarReturnT>
      return_spline(product_degrees);

  // Start calculating the new control points
  for (IndexingType i{}; i < return_spline.GetNumberOfControlPoints(); i++) {
    return_spline.CombineControlPointsForProduct_(
        (*this), rhs, return_spline.LocalToGlobalIndex(i),
        static_cast<ScalarReturnT>(1.));
  }
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                       ScalarRHS>& inner_function) const {
  // Check if the inner function fits into the unit cube
  if (!inner_function.FitsIntoUnitCube()) {
    utils::Logger::TerminatingError(
        "Inner function must fit unit cube to "
        "perform composition");
  }

  // Check Dimensions (Scalar Value splines to be treaded as kSpDim = 1)
  if constexpr (std::is_arithmetic_v<PointTypeRHS>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");
  }
  return ComposeNumeratorSpline(inner_function);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ComposeSensitivity(
        const BezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                           ScalarRHS>& inner_function) const {
  // Check if the inner function fits into the unit cube
  if (!inner_function.FitsIntoUnitCube()) {
    utils::Logger::TerminatingError(
        "Inner function must fit unit cube to "
        "perform composition");
  }

  // Check Dimensions (Scalar Value splines to be treaded as kSpDim = 1)
  if constexpr (std::is_arithmetic_v<PointTypeRHS>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");
  }
  return ComposeNumeratorSensitivity(inner_function);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const RationalBezierSpline<parametric_dimension_inner_spline, PointTypeRHS,
                               ScalarRHS>& inner_function) const {
  // Check if the inner function fits into the unit cube
  if (!inner_function.FitsIntoUnitCube()) {
    utils::Logger::TerminatingError(
        "Inner function must fit unit cube to "
        "perform composition");
  }

  // Check Dimensions (Scalar Value splines to be treaded as kSpDim = 1)
  if constexpr (std::is_arithmetic_v<PointTypeRHS>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");
  }

  return RationalBezierSpline{// Numerator
                              ComposeNumeratorSpline(inner_function),
                              // Denominator
                              ComposeDenominatorSpline(inner_function)};
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ComposeSensitivity(
        const RationalBezierSpline<parametric_dimension_inner_spline,
                                   PointTypeRHS, ScalarRHS>& inner_function)
        const {
  // Type declarations
  using ReturnSplineType =
      RationalBezierSpline<parametric_dimension_inner_spline,
                           decltype(ScalarType{} * ScalarRHS{}),
                           decltype(ScalarType{} * ScalarRHS{})>;

  // Check if the inner function fits into the unit cube
  if (!inner_function.FitsIntoUnitCube()) {
    utils::Logger::TerminatingError(
        "Inner function must fit unit cube to "
        "perform composition");
  }

  // Check Dimensions (Scalar Value splines to be treaded as kSpDim = 1)
  if constexpr (std::is_arithmetic_v<PointTypeRHS>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(PointTypeRHS::kSpatialDimension == parametric_dimension,
                  "Dimension mismatch");
  }

  // Precomputed required values - Numerator
  const auto numerator_splines = ComposeNumeratorSensitivity(inner_function);
  // Precomputed required values - Denominator
  const auto denominator = ComposeDenominatorSpline(inner_function);

  // Init return value
  std::vector<ReturnSplineType> return_value{};
  return_value.reserve(GetNumberOfControlPoints());

  // Create corresponding list of rational splines
  for (IndexingType i{}; i < GetNumberOfControlPoints(); i++) {
    return_value.push_back(ReturnSplineType{numerator_splines[i], denominator});
  }

  return return_value;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename SplineType>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::Compose(
    const std::vector<SplineType>& inner_function_group) const {
  // Ensure Call from BezierType
  static_assert(utils::type_traits::isBezierSpline_v<SplineType> ||
                    utils::type_traits::isRationalBezierSpline_v<SplineType>,
                "Implemented only for (rational) bezier spline-types.");

  /// Start the composition with the current spline
  // Initialize return value
  using ScalarReturnT =
      decltype(ScalarType{} * typename SplineType::ScalarType_{});
  const IndexingType number_of_splines{inner_function_group.size()};
  using ReturnType = std::vector<std::conditional_t<
      utils::type_traits::isRationalBezierSpline_v<SplineType>,
      RationalBezierSpline<SplineType::kParametricDimensions, PhysicalPointType,
                           ScalarReturnT>,
      BezierSpline<SplineType::kParametricDimensions, PhysicalPointType,
                   ScalarReturnT>>>;
  ReturnType composition_group;
  composition_group.reserve(number_of_splines);

  // Check if the inner function fits into the unit cube
  if (!inner_function_group.FitsIntoUnitCube()) {
    utils::Logger::TerminatingError(
        "Inner function must fit unit cube to "
        "perform composition");
  }

  // Check Dimensions
  static_assert(
      SplineType::PhysicalPointType_::kSpatialDimension == parametric_dimension,
      "Dimension mismatch");

  for (IndexingType i_spline{}; i_spline < number_of_splines; i_spline++) {
    composition_group.push_back(Compose(inner_function_group[i_spline]));
  }
  return composition_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <bool compute_sensitivities, typename SplineType>
constexpr auto
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ComposeNumeratorHelper(const SplineType& inner_function) const {
  // Compile only if SplineType is (Rational) Bezier Spline
  static_assert(utils::type_traits::isBezierSpline_v<SplineType> ||
                    utils::type_traits::isRationalBezierSpline_v<SplineType>,
                "ComposeNumerator only compiles for (rational) Beziers.");

  // Deduct Types
  using PointTypeInnerFunction = typename SplineType::PhysicalPointType_;
  using ScalarTypeInnerFunction = typename SplineType::ScalarType_;
  constexpr const std::size_t parametric_dimension_inner_spline =
      SplineType::kParametricDimensions;

  // Sanity check
  if constexpr (std::is_arithmetic_v<PointTypeInnerFunction>) {
    static_assert(parametric_dimension == 1, "Dimension mismatch");
  } else {
    static_assert(
        PointTypeInnerFunction::kSpatialDimension == parametric_dimension,
        "Dimension mismatch");
  }

  // Function is private assuming all checks have already been performed
  /// Start the composition with the current spline
  // Initialize return value
  using ScalarReturnT = decltype(ScalarType_{} * ScalarTypeInnerFunction{});
  const IndexingType sum_of_degrees_outer_spline = std::accumulate(
      degrees.begin(), degrees.end(), static_cast<IndexingType>(0));

  // New degrees
  std::array<IndexingType, parametric_dimension_inner_spline> new_degrees{
      inner_function.GetDegrees()};
  for (IndexingType i{}; i < parametric_dimension_inner_spline; i++)
    new_degrees[i] *= sum_of_degrees_outer_spline;

  // If the sensitivities are to be computed, return a vector with the size of
  // control points and the same degrees, representing the respective
  // derivatives. This falls down to computing the basis function splines of the
  // composition (without the multiplication of the control point values). If
  // the composition is computed directly, the individual splines have top be
  // summed up. Initialize return value
  auto return_value = [&]() {
    if constexpr (compute_sensitivities) {
      using ReturnSplineType = BezierSpline<SplineType::kParametricDimensions,
                                            ScalarReturnT, ScalarReturnT>;
      std::vector<ReturnSplineType> return_value_{};
      // Return one spline for every control point of the outer spline
      return_value_.reserve(GetNumberOfControlPoints());
      return return_value_;
    } else {
      using ReturnSplineType = BezierSpline<SplineType::kParametricDimensions,
                                            PhysicalPointType, ScalarReturnT>;
      ReturnSplineType return_value_{new_degrees};
      return return_value_;
    }
  }();

  // Precompute the "Bernstein-Basis-Splines", similar to Composition of
  // Polynomial splines. This is performed in 2 stages, first the splines are
  // extracted and temporarily stored, then the basis functions are added up
  /*
   * 1. Extract splines, to get the individual physical coordinates
   *    x_i (xi, eta, ...). These are stored temporarily in a local variable
   * 2. Precompute the "Basis function"-splines, i.e., the spline functions
   *    inserted into the definition of the Bernstein polynomials. This might be
   *    storage inefficient, but reduces the total number of multiplications
   *    required.
   */
  std::array<std::vector<BezierSpline<parametric_dimension_inner_spline,
                                      ScalarReturnT, ScalarReturnT>>,
             parametric_dimension>
      bernstein_basis_splines{};
  for (std::size_t i_outer_parametric_dimension{};
       i_outer_parametric_dimension < parametric_dimension;
       i_outer_parametric_dimension++) {
    // Store current parametric dimensions entry as a reference
    auto& bernstein_spline =
        bernstein_basis_splines[i_outer_parametric_dimension];
    bernstein_spline.resize(degrees[i_outer_parametric_dimension] + 1);

    // Extract the current dimension of the spline and its inverse
    auto extracted_dimension_multiple =
        BezierSpline<parametric_dimension_inner_spline, ScalarReturnT,
                     ScalarReturnT>{
            std::array<std::size_t,
                       parametric_dimension_inner_spline>{},  // initializes
                                                              // to zero
            std::vector<ScalarReturnT>{static_cast<ScalarReturnT>(1.)}};
    auto inverse_extracted_dimension_mulitple = extracted_dimension_multiple;

    const auto extracted_dimension =
        inner_function.ExtractDimension(i_outer_parametric_dimension);

    // If the inner function is a rational spline, the weight is multiplied into
    // the first paranthesis to performe the composition. If, however, the inner
    // function is a polynomial spline, the weight function is just 1. In order
    // to keep the variable constant, we use a lamda expression for assignment
    // here
    // For a polynomial inner funciton c(t)
    // $$ \sum (1 - c(t))^i c(t)^{n-i} $$
    // else for rational splines (c(t) / w(t)) the numerator reads
    // $$ \sum (w(t) - c(t)) ^i (c(t)) ^{n-i}$$
    const auto inverse_extracted_dimension = [&]() {
      if constexpr (utils::type_traits::isRationalBezierSpline_v<SplineType>) {
        return inner_function.GetWeightFunctionSpline() - extracted_dimension;

      } else {
        return static_cast<ScalarTypeInnerFunction>(1.) - extracted_dimension;
      }
    }();

    // Loop over all possible Bernstein basis functions and assign the n-th
    // entry the extracted dimension to the power n (see "Free Form Surface
    // Analysis Using a Hybrid of Symbolic and Numeric Computation" - Elber -
    // Eq. 6.2)
    // 1. Multiply c(t) ^ i - terms
    bernstein_spline[0] = extracted_dimension_multiple;
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim <= degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      extracted_dimension_multiple =
          extracted_dimension_multiple * extracted_dimension;
      bernstein_spline[i_ctps_in_par_dim] = extracted_dimension_multiple;
    }
    // 2. Multiply (1 - c(t)) ^ (n-i) - terms
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim <= degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      inverse_extracted_dimension_mulitple =
          inverse_extracted_dimension_mulitple * inverse_extracted_dimension;
      bernstein_spline[degrees[i_outer_parametric_dimension] -
                       i_ctps_in_par_dim] =
          bernstein_spline[degrees[i_outer_parametric_dimension] -
                           i_ctps_in_par_dim] *
          inverse_extracted_dimension_mulitple;
    }
    // 3. Multiply coefficients
    for (std::size_t i_ctps_in_par_dim{1};
         i_ctps_in_par_dim < degrees[i_outer_parametric_dimension];
         i_ctps_in_par_dim++) {
      bernstein_spline[degrees[i_outer_parametric_dimension] -
                       i_ctps_in_par_dim] *=
          utils::FastBinomialCoefficient::choose(
              degrees[i_outer_parametric_dimension], i_ctps_in_par_dim);
    }
  }

  /*
   * @todo The current loop can be made more efficient by precomputing spline
   * products
   *
   * for the three-dimensional case C = T_{ijk} Bi Bj Bk
   * Looping over all i j k (as is currently done) is very inefficient, as many
   * splines Products are computed multiple times (ijk). If the product can be
   * organized as follows, the computational cost can significantly be reduced.
   * for i:
   *    for j:
   *        Bij = Bi * Bj
   *        for k:
   *            Bijk = Bij * Bk
   *
   * In this example for every k -- (k-1) fewer spline products are computed as
   * Bij can be reused for every iteration of k. This effect becomes more
   * important for parametric dimensions of the deformation function
   */
  // Loop over control points of the outer spline to start the interpolation
  for (unsigned int i_control_point_outer_spline{};
       i_control_point_outer_spline < GetNumberOfControlPoints();
       i_control_point_outer_spline++) {
    // Retrieve indices
    const auto outer_spline_ctps_indices =
        LocalToGlobalIndex(i_control_point_outer_spline);

    // Start interpolation
    if constexpr (compute_sensitivities) {
      return_value.push_back(
          bernstein_basis_splines[0][outer_spline_ctps_indices[0]]);
      for (std::size_t i_outer_parametric_dimension{1};
           i_outer_parametric_dimension < parametric_dimension;
           i_outer_parametric_dimension++) {
        return_value.back() =
            return_value.back() *
            bernstein_basis_splines
                [i_outer_parametric_dimension]
                [outer_spline_ctps_indices[i_outer_parametric_dimension]];
      }
    } else {
      // Initialize current factor
      auto bernstein_basis_tensor_product =
          bernstein_basis_splines[0][outer_spline_ctps_indices[0]];
      for (std::size_t i_outer_parametric_dimension{1};
           i_outer_parametric_dimension < parametric_dimension;
           i_outer_parametric_dimension++) {
        bernstein_basis_tensor_product =
            bernstein_basis_tensor_product *
            bernstein_basis_splines
                [i_outer_parametric_dimension]
                [outer_spline_ctps_indices[i_outer_parametric_dimension]];
      }

      // Check if the bernstein_basis_tensor_product has the degrees (and number
      // of control_points) as expected
      assert(new_degrees == bernstein_basis_tensor_product.GetDegrees());

      // Now that all factors have been precomputed we can determine the
      // control-points of the resulting spline
      for (IndexingType i_control_point_composition{};
           i_control_point_composition <
           return_value.GetNumberOfControlPoints();
           i_control_point_composition++) {
        return_value.control_points[i_control_point_composition] +=
            bernstein_basis_tensor_product
                .control_points[i_control_point_composition] *
            control_points[i_control_point_outer_spline];
      }
    }
  }

  return return_value;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <std::size_t parametric_dimension_inner_spline, typename PointTypeRHS,
          typename ScalarRHS>
constexpr BezierSpline<parametric_dimension_inner_spline,
                       decltype(ScalarType{} * ScalarRHS{}),
                       decltype(ScalarType{} * ScalarRHS{})>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ComposeDenominatorSpline(
        const RationalBezierSpline<parametric_dimension_inner_spline,
                                   PointTypeRHS, ScalarRHS>& inner_function)
        const {
  // Calculate cumulated sum of outer degrees
  const IndexingType sum_of_degrees_outer_spline = std::accumulate(
      degrees.begin(), degrees.end(), static_cast<IndexingType>(0));
  return (inner_function.GetWeightFunctionSpline())
      .RaisePower(sum_of_degrees_outer_spline);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::vector<
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>
BezierSpline<parametric_dimension, PhysicalPointType,
             ScalarType>::SplitAtPosition(const ScalarType& splitting_plane,
                                          const IndexingType
                                              splitting_dimension) const {
  utils::Logger::ExtendedInformation(
      "Splitting spline along parametric dimension " +
      std::to_string(splitting_dimension) + " at position " +
      std::to_string(splitting_plane));
  /*
   * Notation:
   * A group the number of all starting points that are subsequent, i.e. the
   * size of a group is the number of subsequent 0s. E.g. in the example:
   * Group size of index 0 is 1, group size of index 1 is 3, etc.. Group
   * offset is thus the distance between two group-starting-points. As before,
   * the index offset is the offset to the next layer with the same
   * corresponding local indices in the other parametric dimensions.
   *
   *    local indices : global index
   *    0   0   0 : 0   -   -
   *    1   0   0 : 1   |   |  group 0 for index 1
   *    2   0   0 : 2   |   -
   *    0   1   0 : 3   |
   *    1   1   0 : 4   | group offset
   *    2   1   0 : 5   |
   *    0   0   1 : 6   -   -
   *    1   0   1 : 7       |  group 1 for index 1
   *    2   0   1 : 8       -
   *    0   1   1 : 9
   *    1   1   1 : 10
   *    2   1   1 : 11
   */

  // Check bounds (as default knot-vector is [0,0,..,1,1,...])
  if (!(splitting_plane > 0. && splitting_plane < 1.)) {
    utils::Logger::TerminatingError("Split Plane is at " +
                                    std::to_string(splitting_plane) +
                                    ", which is out of bounds (0,1)");
  }

  // Initialize return values
  BezierSpline first_part{(*this)}, second_part{(*this)};
  using ReturnType = std::vector<
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>;
  ReturnType spline_group_split{first_part, second_part};

  // References for simplified access
  auto& ctps_first = spline_group_split[0].control_points;
  auto& ctps_second = spline_group_split[1].control_points;

  // Control points per layer normal to parametric direction
  IndexingType ctps_per_layer =
      GetNumberOfControlPoints() / (degrees[splitting_dimension] + 1);

  // Evaluation factors
  const ScalarType_ factor = splitting_plane;
  const ScalarType_ inv_factor = 1. - splitting_plane;

  // Auxiliary index offsets
  const IndexingType starting_points_per_index_group =
      index_offsets[splitting_dimension];
  const IndexingType group_offset =
      index_offsets[splitting_dimension] * (degrees[splitting_dimension] + 1);
  const IndexingType index_offset = index_offsets[splitting_dimension];

  // Loop over all starting points
  for (IndexingType start_point_id{0}; start_point_id < ctps_per_layer;
       start_point_id++) {
    // Locate position of starting point
    const IndexingType group_id =
        start_point_id / starting_points_per_index_group;
    const IndexingType id_in_group =
        start_point_id % starting_points_per_index_group;
    const IndexingType global_start_id = group_id * group_offset + id_in_group;

    // Going forward on first part and taking last value to second part
    for (IndexingType i{1}; i <= degrees[splitting_dimension]; i++) {
      // Reversed order to prohibit data overwrite
      for (IndexingType j{degrees[splitting_dimension]}; j >= i; j--) {
        ctps_first[global_start_id + j * index_offset] =
            ctps_first[global_start_id + (j - 1) * index_offset] * inv_factor +
            ctps_first[global_start_id + j * index_offset] * factor;
      }
      ctps_second[global_start_id +
                  (degrees[splitting_dimension] - i) * index_offset] =
          ctps_first[global_start_id +
                     degrees[splitting_dimension] * index_offset];
    }
  }
  return spline_group_split;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::vector<
    BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>
BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    SplitAtPosition(const std::vector<ScalarType>& splitting_planes,
                    const IndexingType splitting_dimension) const {
  // Init return value
  using ReturnType = std::vector<BezierSpline>;
  ReturnType return_value;
  return_value.reserve(splitting_planes.size() + 1);

  // Check if the split plane vector has increasing values
  if (!std::is_sorted(splitting_planes.begin(), splitting_planes.end())) {
    utils::Logger::TerminatingError("Split Plane values must be increasing");
  }

  // If vector is empty return empty group containing only the Spline itself
  if (splitting_planes.size() == 0) {
    return_value.push_back(*this);
    return return_value;
  }

  // Loop over all split planes and assign the first part to the vector,
  // then split the remaining part at the relative required position
  auto split_spline = SplitAtPosition(splitting_planes[0], splitting_dimension);
  return_value.push_back(split_spline[0]);
  for (std::size_t i_split{1}; i_split < splitting_planes.size(); i_split++) {
    // Calculate the position of the split position in the remaining spline
    const double split_in_remaining_part =
        (splitting_planes[i_split] - splitting_planes[i_split - 1]) /
        (1. - splitting_planes[i_split - 1]);
    // Split remainder
    split_spline = split_spline[1].SplitAtPosition(split_in_remaining_part,
                                                   splitting_dimension);
    return_value.push_back(split_spline[0]);
  }
  // Assign last section to last element
  return_value.push_back(split_spline[1]);
  return return_value;
}
