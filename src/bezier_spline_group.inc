#ifndef SRC_BEZIER_SPLINE_GROUP_INC
#define SRC_BEZIER_SPLINE_GROUP_INC

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr bool BezierSplineGroup<parametric_dimension, PhysicalPointType,
                                 ScalarType>::fits_unit_cube() const {
  for (auto spline : (*this)) {
    if (!spline.fits_unit_cube()) {
      return false;
    }
  }
  return true;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType BezierSplineGroup<
    parametric_dimension, PhysicalPointType, ScalarType>::maximum() const {
  PhysicalPointType maximum_corner = (*this)[0].maximum();
  for (IndexingType i{1}; i < this->size(); i++) {
    const auto maximum_corner_i = (*this)[i].maximum();
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      maximum_corner[i_dim] =
          std::max(maximum_corner[i_dim], maximum_corner_i[i_dim]);
    }
  }
  return maximum_corner;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr PhysicalPointType BezierSplineGroup<
    parametric_dimension, PhysicalPointType, ScalarType>::minimum() const {
  PhysicalPointType minimum_corner = (*this)[0].minimum();
  for (IndexingType i{1}; i < this->size(); i++) {
    const auto minimum_corner_i = (*this)[i].minimum();
    for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
         i_dim++) {
      minimum_corner[i_dim] =
          std::min(minimum_corner[i_dim], minimum_corner_i[i_dim]);
    }
  }
  return minimum_corner;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>&
BezierSplineGroup<parametric_dimension, PhysicalPointType,
                  ScalarType>::fit_to_unit_cube() {
  const auto maximum_corner = maximum();
  const auto minimum_corner = minimum();
  auto scale_vector = maximum_corner - minimum_corner;
  for (IndexingType i_dim{}; i_dim < PhysicalPointType::kSpatialDimension;
       i_dim++) {
    scale_vector[i_dim] = static_cast<ScalarType>(1) / scale_vector[i_dim];
  }
  for (IndexingType i{1}; i < this->size(); i++) {
    (*this)[i].transpose_and_scale(-minimum_corner, scale_vector);
  }
  return (*this);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::compose(
    const SplineBaseType& inner_function) const {
  // Initiate return value
  const auto number_of_splines_in_group = this->size();
  BezierSplineGroup composed_group{number_of_splines_in_group};

  // Compose each part of the structure individually
  for (IndexingType i{}; i < number_of_splines_in_group; i++) {
    composed_group[i] = (*this)[i].compose(inner_function);
  }
  return composed_group;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
constexpr BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>::compose(
    const BezierSplineGroup& inner_function_group) const {
  // Initiate return value
  const auto number_of_splines_in_outer_group = this->size();
  const auto number_of_splines_in_inner_group = inner_function_group.size();

  BezierSplineGroup composed_group{number_of_splines_in_outer_group *
                                   number_of_splines_in_inner_group};

  // Compose each part of the structure individually
  for (IndexingType i{}; i < number_of_splines_in_inner_group; i++) {
    const auto i_composition = compose(inner_function_group[i]);
    for (IndexingType j{}; j < number_of_splines_in_outer_group; j++) {
      composed_group[i * number_of_splines_in_outer_group + j] =
          i_composition[j];
    }
  }
  return composed_group;
}

#endif  // SRC_BEZIER_SPLINE_GROUP_INC
