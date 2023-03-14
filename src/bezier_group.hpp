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

#ifndef SRC_BEZIER_GROUP_HPP
#define SRC_BEZIER_GROUP_HPP

#include <cassert>
#include <vector>

#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/rational_bezier_spline.hpp"
#include "bezman/src/utils/logger.hpp"
#include "bezman/src/utils/type_traits/is_bezier_spline.hpp"
#include "bezman/src/utils/type_traits/is_rational_bezier_spline.hpp"

namespace bezman {

/*
 * Group of Bezier splines with same parametric dimension and PointType
 *
 * Predefining the PointType and the parametric dimension is somewhat
 * restricting, but avoids overcomplicated CRTP implementation. Using tuples
 * woud restrict the use of this function to compile time descriptions and
 * prevents use of export import routines.
 *
 * Class inherits from std::vector. This facilitates the use of its methods in
 * later applications. No operator[] overload required etc.
 *
 * In the microstructure use case, BezierGroup can be used both ways,
 * either as the Microtile, or as the deformation function patches, that form
 * the
 *
 * @tparam SplineType : BezierSpline or RationalBezierSpline
 *
 * All other SplineTypes will fail
 */
template <typename SplineType>
class BezierGroup : public std::vector<SplineType> {
  // Make sure that this class is only build using (Rational) Bezier Spline
  // Types
  static_assert(utils::type_traits::isBezierSpline_v<SplineType> ||
                    utils::type_traits::isRationalBezierSpline_v<SplineType>,
                "BezierGroup may only contain (rational) bezier splines.");

 public:
  using SplineType_ = SplineType;
  using PhysicalPointType_ = typename SplineType::PhysicalPointType_;
  using ScalarType_ = typename SplineType::ScalarType_;

  template <typename SplineTypeRHS>
  using ComposedType = decltype(SplineType{}.Compose(SplineTypeRHS{}));
  using ScalarSplineType = decltype(SplineType{}.ExtractDimension(0));
  //   std::conditional_t <
  //   utils::type_traits::isRationalBezierSpline_v<SplineType> ||
  //       utils::type_traits::isRationalBezierSpline_v<SplineTypeRHS>,
  //       RationalBezierSpline < SplineTypeRHS::kParametricDimensions, typename
  //       SplineType

 private:
  using IndexingType = std::size_t;
  using BaseVector = std::vector<SplineType>;

 public:
  /// Parametric dimension of SplineType
  constexpr static IndexingType kParametricDimensions =
      SplineType::kParametricDimensions;

  /// Default constructor (to profit from std::vectors implementations)
  constexpr BezierGroup() = default;

  /// Default constructor (to profit from std::vectors implementations)
  template <typename IntegralType, typename = typename std::enable_if_t<
                                       std::is_integral_v<IntegralType>>>
  constexpr BezierGroup(const IntegralType &init_size)
      : std::vector<SplineType_>{init_size} {}

  /// Copy constructor
  constexpr BezierGroup(const BezierGroup &) = default;

  /// Initializer list overload
  template <typename... Splines>
  constexpr BezierGroup(const Splines &...splines)
      : BaseVector{static_cast<SplineType>(splines)...} {}

  /// Check if group fits unit cube
  constexpr bool FitsIntoUnitCube() const;

  /// Maximum
  constexpr PhysicalPointType_ MaximumCorner() const;

  /// Minimum corner
  constexpr PhysicalPointType_ MinimumCorner() const;

  /// Fit to unit_cube
  constexpr BezierGroup &FitIntoUnitCube();

  /// Compose with single Spline
  template <typename SplineTypeRHS>
  constexpr auto Compose(const SplineTypeRHS &inner_function) const;

  /// Compose with Splinegroup
  template <typename SplineTypeRHS>
  constexpr auto Compose(
      const BezierGroup<SplineTypeRHS> &inner_function_group) const;

  /// Add two Bezier Spline Groups Component wise to the current Group
  constexpr BezierGroup &AddComponentwise(const BezierGroup &rhs);

  /// Add two Bezier Spline Groups Component wise to the current Group
  template <typename SplineTypeRHS>
  constexpr auto MultiplyComponentwise(
      const BezierGroup<SplineTypeRHS> &rhs) const;

  /// Calculate the derivative of all components and return in a new group
  constexpr BezierGroup DerivativeWRTParametricDimension(
      const IndexingType par_dim) const;

  /// Extract Dimensions Component-wise
  constexpr BezierGroup<ScalarSplineType> ExtractDimension(
      const IndexingType &par_dim) const;

  /// + Operator for concatenation
  constexpr BezierGroup operator+(const BezierGroup &rhs) const;

  /// + Operator for concatenation
  constexpr BezierGroup &operator+=(const BezierGroup &rhs);

  /// + Operator for translation
  constexpr BezierGroup operator+(const PhysicalPointType_ &translation) const;

  /// + Operator for translation
  constexpr BezierGroup &operator+=(const PhysicalPointType_ &translation);
};  // namespace bezman

#include "bezman/src/bezier_group.inc"

}  // namespace bezman

#endif  // SRC_BEZIER_GROUP_HPP
