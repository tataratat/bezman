#ifndef UTILS_TYPE_TRAITS_IS_POINT_HPP
#define UTILS_TYPE_TRAITS_IS_POINT_HPP

#include "bezierManipulation/src/point.hpp"

namespace beziermanipulation::utils::type_traits
{
  
/// Checker if a template type is an instance of Point
template <typename PointType>
struct isPoint {
  constexpr static bool value = false;
};

/// Checker if a template type is an instance of Point
template <std::size_t spatial_dimension, typename BaseType>
struct isPoint<Point<spatial_dimension, BaseType>> {
  constexpr static bool value = true;
};

/// Alias std naming conform
template <typename PointType>
inline constexpr bool isPoint_v = isPoint<PointType>::value;

} // namespace beziermanipulation::utils::type_traits


#endif  // UTILS_TYPE_TRAITS_IS_POINT_HPP
