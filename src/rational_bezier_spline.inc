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
template <typename... T>
constexpr PhysicalPointType
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::ControlPoint(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  const std::array<IndexingType, parametric_dimension> indices{
      (static_cast<IndexingType>(index), ...)};
  return ControlPoint(indices);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    ControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) const {
  return weighted_spline_.ControlPoint(index) /
         weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::WeightedControlPoint(const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weighted_spline_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::WeightedControlPoint(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weighted_spline_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    WeightedControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) {
  return weighted_spline_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const PhysicalPointType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    WeightedControlPoint(
        const std::array<IndexingType, parametric_dimension>& index) const {
  return weighted_spline_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::Weight(const T... index) {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weight_function_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename... T>
constexpr const ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::Weight(const T... index) const {
  static_assert(sizeof...(T) == parametric_dimension,
                "Unspecified number of indices.");
  return weight_function_.ControlPoint(index...);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    Weight(const std::array<IndexingType, parametric_dimension>& index) {
  return weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr const ScalarType&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    Weight(const std::array<IndexingType, parametric_dimension>& index) const {
  return weight_function_.ControlPoint(index);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType RationalBezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::Evaluate(const PointTypeParametric_& par_coords) const {
  return weighted_spline_.Evaluate(par_coords) /
         weight_function_.Evaluate(par_coords);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    BasisFunctions(const PointTypeParametric_& par_coords) const {
  // Check sanity
  assert(CheckSplineCompatibility());
  // Evaluiere Basis Functions
  std::array<std::vector<ScalarType>, parametric_dimension> factors =
      weighted_spline_.BasisFunctions(par_coords);

  // Initialize values
  ScalarType weight_eval{};
  weight_function_.AddUpContributionsToControlPointVector_(
      weight_eval, factors, static_cast<ScalarType>(1.));

  // Add up the individual contributions
  const ScalarType inv_weight_function =
      static_cast<ScalarType>(1.) / weight_eval;
  for (IndexingType par_dim{0}; par_dim < parametric_dimension; par_dim++) {
    for (IndexingType i_basis_function{0};
         i_basis_function < factors[par_dim].size(); i_basis_function++) {
      factors[par_dim][i_basis_function] *= inv_weight_function;
    }
  }
  return factors;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr std::array<std::vector<ScalarType>, parametric_dimension>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    PolynomialBasisFunctions(const PointTypeParametric_& par_coords) const {
  // Pass to member
  return weighted_spline_.BasisFunctions(par_coords);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType RationalBezierSpline<
    parametric_dimension, PhysicalPointType,
    ScalarType>::ForwardEvaluate(const PointTypeParametric_& par_coords) const {
  // Check sanity
  assert(CheckSplineCompatibility());

  // Evaluiere Basis Functions
  const auto factors = weighted_spline_.BasisFunctions(par_coords);

  // Initialize values
  PhysicalPointType evaluation_point{};
  ScalarType weight_eval{};

  // Add up the individual contributions
  return weighted_spline_.AddUpContributionsToControlPointVector_(
             evaluation_point, factors, static_cast<ScalarType>(1.)) /
         weight_function_.AddUpContributionsToControlPointVector_(
             weight_eval, factors, static_cast<ScalarType>(1.));
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    EvaluateDerivative(
        const PointTypeParametric_& par_coords,
        const std::array<std::size_t, parametric_dimension>& nth_derivs) const {
  // Define Lambdas for indexing
  // From a single index notation to an array of indices
  auto local_ids_ = [&nth_derivs](const std::size_t req_id)
      -> std::array<std::size_t, parametric_dimension> {
    std::size_t id{req_id};
    std::array<std::size_t, parametric_dimension> local_ids{};
    for (std::size_t i_pd{}; i_pd < parametric_dimension; i_pd++) {
      if (nth_derivs[i_pd] == 0) continue;
      local_ids[i_pd] = id % (nth_derivs[i_pd] + 1);
      id -= local_ids[i_pd];
      id /= nth_derivs[i_pd] + 1;
    }
    return local_ids;
  };
  // From array of indices to single index
  auto global_ids_ =
      [&nth_derivs](
          const std::array<std::size_t, parametric_dimension>& req_derivs)
      -> std::size_t {
    std::size_t id{};
    std::size_t offset{1};
    for (std::size_t i_pd{}; i_pd < parametric_dimension; i_pd++) {
      // assert(req_derivs[i_pd] <= nth_derivs[i_pd]);
      if (nth_derivs[i_pd] == 0) continue;
      id += offset * (req_derivs[i_pd]);
      offset = nth_derivs[i_pd] > 0 ? offset * (nth_derivs[i_pd] + 1) : offset;
    }
    return id;
  };

  // Initialize return type
  const std::size_t number_of_derivs{global_ids_(nth_derivs) + 1};
  // Please remember that the first derivative is not used
  std::vector<PhysicalPointType> derivatives(number_of_derivs);
  std::vector<PhysicalPointType> A_derivatives(number_of_derivs);
  std::vector<ScalarType> w_derivatives(number_of_derivs);

  // Fill all polynomial spline derivatives (and values for id=0)
  for (std::size_t i_deriv{}; i_deriv < number_of_derivs; i_deriv++) {
    const auto req_derivs = local_ids_(i_deriv);
    A_derivatives[i_deriv] =
        weighted_spline_.EvaluateDerivative(par_coords, req_derivs);
    w_derivatives[i_deriv] =
        weight_function_.EvaluateDerivative(par_coords, req_derivs);
  }

  // Precompute inverse of weighted function
  const ScalarType inv_w_fact = static_cast<ScalarType>(1.) / w_derivatives[0];

  // Loop over all lower-order derivatives and assign derivatives-vector
  // Notation follows "The NURBS book" eq. 4.20 (extended for n-d splines)
  for (std::size_t i_deriv{0}; i_deriv < number_of_derivs; i_deriv++) {
    // Retrieve index-wise order of the derivative for current ID
    const auto derivative_order_indexwise_LHS = local_ids_(i_deriv);
    // Assign derivative of Numerator-function
    derivatives[i_deriv] = A_derivatives[i_deriv];
    // Substract all weighted lower-order functions
    for (std::size_t j_deriv{0}; j_deriv < i_deriv; j_deriv++) {
      // Retrieve order of current index
      const auto derivative_order_indexwise_RHS = local_ids_(j_deriv);
      // Precompute Product of binomial coefficients
      std::size_t binom_fact{1};
      for (std::size_t i_pd{}; i_pd < parametric_dimension; i_pd++) {
        binom_fact *= utils::FastBinomialCoefficient::choose(
            derivative_order_indexwise_LHS[i_pd],
            derivative_order_indexwise_RHS[i_pd]);
      }
      // Substract low-order function
      derivatives[i_deriv] -=
          binom_fact * w_derivatives[i_deriv - j_deriv] * derivatives[j_deriv];
    }
    // Finalize
    derivatives[i_deriv] *= inv_w_fact;
  }
  // Return last value
  return derivatives[number_of_derivs - 1];
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    OrderElevateAlongParametricDimension(const IndexingType par_dim) {
  weighted_spline_.OrderElevateAlongParametricDimension(par_dim);
  weight_function_.OrderElevateAlongParametricDimension(par_dim);
  // Check sanity
  assert(CheckSplineCompatibility());
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierGroup<
    RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    SplitAtPosition(const ScalarType& splitting_planes,
                    const IndexingType splitting_dimension) const {
  const auto weighted_splines =
      weighted_spline_.SplitAtPosition(splitting_planes, splitting_dimension);
  const auto weight_functions =
      weight_function_.SplitAtPosition(splitting_planes, splitting_dimension);
  return BezierGroup<RationalBezierSpline>{
      RationalBezierSpline{weighted_splines[0], weight_functions[0]},
      RationalBezierSpline{weighted_splines[1], weight_functions[1]}};
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierGroup<
    RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    SplitAtPosition(const std::vector<ScalarType>& splitting_planes,
                    const IndexingType splitting_dimension) const {
  // Init return value
  using ReturnType = BezierGroup<RationalBezierSpline>;
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

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
    DerivativeWRTParametricDimension(const IndexingType par_dim) const {
  // Implemented using Quotient-rule
  // Note that the numerator spline needs to be order elevated to match
  // denominator's order
  return RationalBezierSpline(
      // u' * w - w' * u
      (weighted_spline_.DerivativeWRTParametricDimension(par_dim) *
           weight_function_ -
       weight_function_.DerivativeWRTParametricDimension(par_dim) *
           weighted_spline_)
          .OrderElevateAlongParametricDimension(par_dim),
      // v^2
      weight_function_ * weight_function_);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator*=(const ScalarType& scalar) {
  // Multiplication only reflects in the weighted spline
  weighted_spline_ *= scalar;
  // Check sanity
  assert(CheckSplineCompatibility());
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator*(const ScalarType& scalar) const {
  // Multiplication only reflects in the weighted spline
  return RationalBezierSpline(weighted_spline_ * scalar, weight_function_);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+=(const PhysicalPointType& point_shift) {
  // Add up contributions to control points with weights
  for (IndexingType i_ctps{};
       i_ctps < weighted_spline_.GetNumberOfControlPoints(); i_ctps++) {
    weighted_spline_.control_points[i_ctps] +=
        point_shift * weight_function_.control_points[i_ctps];
  }
  // Check sanity
  assert(CheckSplineCompatibility());
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+(const PhysicalPointType& point_shift) const {
  // Pass to operator +=
  RationalBezierSpline return_spline{weighted_spline_, weight_function_};
  return_spline += point_shift;
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-=(const PhysicalPointType& point_shift) {
  // Add up contributions to control points with weights
  for (IndexingType i_ctps{};
       i_ctps < weighted_spline_.GetNumberOfControlPoints(); i_ctps++) {
    weighted_spline_.control_points[i_ctps] -=
        point_shift * weight_function_.control_points[i_ctps];
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-(const PhysicalPointType& point_shift) const {
  // Pass to operator +=
  RationalBezierSpline return_spline{weighted_spline_, weight_function_};
  return_spline -= point_shift;
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::operator-() const {
  // Pass to operator +=
  RationalBezierSpline return_spline{-weighted_spline_, weight_function_};
  return return_spline;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} + PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  if (weight_function_ == rhs.weight_function_) {
    return RationalBezierSpline{weighted_spline_ + rhs.weighted_spline_,
                                weight_function_};
  } else {
    return RationalBezierSpline{weighted_spline_ * rhs.weight_function_ +
                                    weight_function_ * rhs.weighted_spline_,
                                weight_function_ * rhs.weight_function_};
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator+=(const RationalBezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& rhs) {
  if (weight_function_ == rhs.weight_function_) {
    weighted_spline_ += rhs.weighted_spline_;
    return (*this);
  } else {
    weighted_spline_ = weighted_spline_ * rhs.weight_function_;
    weighted_spline_ += weight_function_ * rhs.weighted_spline_;
    weight_function_ = weight_function_ * rhs.weight_function_;
    // Check sanity
    assert(CheckSplineCompatibility());
    return (*this);
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} - PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  if (weight_function_ == rhs.weight_function_) {
    return RationalBezierSpline{weighted_spline_ - rhs.weighted_spline_,
                                weight_function_};
  } else {
    return RationalBezierSpline{weighted_spline_ * rhs.weight_function_ -
                                    weight_function_ * rhs.weighted_spline_,
                                weight_function_ * rhs.weight_function_};
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr RationalBezierSpline<parametric_dimension, PhysicalPointType,
                               ScalarType>&
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator-=(const RationalBezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& rhs) {
  if (weight_function_ == rhs.weight_function_) {
    weighted_spline_ -= rhs.weighted_spline_;
    return (*this);
  } else {
    weighted_spline_ = weighted_spline_ * rhs.weight_function_;
    weighted_spline_ -= weight_function_ * rhs.weighted_spline_;
    weight_function_ = weight_function_ * rhs.weight_function_;
    // Check sanity
    assert(CheckSplineCompatibility());
    return (*this);
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
template <typename PointTypeRHS, typename ScalarRHS>
constexpr RationalBezierSpline<parametric_dimension,
                               decltype(PhysicalPointType{} * PointTypeRHS{}),
                               decltype(ScalarType{} * ScalarRHS{})>
RationalBezierSpline<parametric_dimension, PhysicalPointType, ScalarType>::
operator*(const RationalBezierSpline<parametric_dimension, PointTypeRHS,
                                     ScalarRHS>& rhs) const {
  return RationalBezierSpline<parametric_dimension,
                              decltype(PhysicalPointType{} * PointTypeRHS{}),
                              decltype(ScalarType{} * ScalarRHS{})>{
      weighted_spline_ * rhs.weighted_spline_,
      weight_function_ * rhs.weight_function_};
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSpline<parametric_dimension, ScalarType, ScalarType>
RationalBezierSpline<parametric_dimension, PhysicalPointType,
                     ScalarType>::ExtractDimension(const std::size_t& dimension)
    const {
  return weighted_spline_.ExtractDimension(dimension);
}
