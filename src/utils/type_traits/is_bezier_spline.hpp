#ifndef UTILS_TYPE_TRAITS_IS_BEZIER_SPLINE_HPP
#define UTILS_TYPE_TRAITS_IS_BEZIER_SPLINE_HPP

#include "bezman/src/bezier_spline.hpp"

namespace bezman::utils::type_traits {

/// Checker if a template type is an instance of BezierSpline
template <typename SplineType>
struct isBezierSpline {
  constexpr static bool value = false;
};

/// Checker if a template type is an instance of BezierSpline
template <std::size_t parametric_dimension, typename PhysicalSplineType,
          typename ScalarType>
struct isBezierSpline<BezierSpline<spatial_dimension, BaseType>> {
  constexpr static bool value = true;
};

/// Alias std naming conform
template <typename SplineType>
constexpr bool isBezierSpline_v = isBezierSpline<SplineType>::value;

}  // namespace bezman::utils::type_traits

#endif  // UTILS_TYPE_TRAITS_IS_BEZIER_SPLINE_HPP
