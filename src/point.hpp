#ifndef SRC_POINT_HPP
#define SRC_POINT_HPP

#include <array>
#include <iomanip>
#include <iostream>

namespace beziermanipulation {

// MINIMALISTIC COPY FROM CAMPIGA POINTS
/*
 * NEEDS DOCUMENTATION @tbd
 */
template <unsigned int space_dim, typename BaseType = double>
class Point : public std::array<BaseType, space_dim> {
 public:
   constexpr static unsigned int kSpatialDimension = space_dim;


   
   constexpr Point(const Point &) = default;

   constexpr Point() = default;

   template <typename... scalar>
   explicit constexpr Point(const scalar &...coords)
       : std::array<BaseType, space_dim>{coords...}
   {
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

  template <typename Scalar>
  decltype(Scalar{} * BaseType{}) operator*(
      const Point<space_dim, Scalar>& point) const {
    using ScalarReturn = decltype(Scalar{} * BaseType{});
    ScalarReturn result{};
    for (unsigned int i{}; i < space_dim; i++) {
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
};  // namespace std::array<BaseType,space_dim>
}  // namespace beziermanipulation

#endif  // SRC_POINT_HPP
