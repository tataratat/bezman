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