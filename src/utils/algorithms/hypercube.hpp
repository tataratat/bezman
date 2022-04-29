#ifndef UTILS_TYPE_TRAITS_HYPERCUBE_HPP
#define UTILS_TYPE_TRAITS_HYPERCUBE_HPP

#include <array>

#include "bezierManipulation/src/utils/algorithms/int_power.hpp"

namespace beziermanipulation::utils::algorithms {

/**
 * @brief Provides functions and auxiliary values for hypercubes
 *
 * This can be used to retrieve indices for values in the parametric space
 *
 * @tparam dimension
 */
template <std::size_t dimension>
class HyperCube {
  // Current implementation limited to 2 and 3 dimensions
  static_assert((dimension == 3 || dimension == 2),
                "High-Dimensional and Line Patches not supported");

 public:
  /**
   * @brief Get the Opposite Faces
   *
   * @return constexpr std::array<std::size_t, dimension * 2>
   *                      Array with corresponding faces
   */
  static constexpr std::array<std::size_t, dimension * 2> GetOppositeFaces() {
    if constexpr (dimension == 2) {
      return std::array<std::size_t, dimension * 2>{2, 3, 0, 1};
    } else if constexpr (dimension == 3) {
      // Must be 3D
      return std::array<std::size_t, dimension * 2>{5, 3, 4, 1, 2, 0};
    } else {
      // Should never happen
      return std::array<std::size_t, dimension * 2>{};
    }
  }
  /**
   * @brief Control point indices of associated faces
   *
   * @param degrees
   * @return auto
   */
  static auto ControlPointIndicesOnFace(
      const std::array<std::size_t, dimension>& degrees) {
    static_assert((dimension == 3 || dimension == 2),
                  "High-Dimensional and Line Patches not supported");
    if constexpr (dimension == 2) {
      /* returns:
       *
       *  3 --(2)>- 2
       *  ^         ^
       * (3)   0   (1)
       *  |         |
       *  0 --(0)>- 1
       */
      return std::array<std::array<std::size_t, 2>, dimension * 2>{
          0,
          degrees[0],
          degrees[0],
          (degrees[0] + 1) * (degrees[1] + 1) - 1,
          (degrees[0] + 1) * degrees[1],
          (degrees[0] + 1) * (degrees[1] + 1) - 1,
          0,
          (degrees[0] + 1) * degrees[1]};
    } else {
      // Must be 3D (@todo implement)
      static_assert(dimension == 2, "Not Implemented");
      return std::array<std::array<std::size_t, 4>, dimension * 2>{};
    }
  }
  /**
   * @brief Vertex indices (of corners) on faces
   *
   * @return constexpr std::array<std::size_t, 2^(dim-1)> vertices on faces
   */
  static constexpr std::array<
      std::array<std::size_t, beziermanipulation::utils::algorithms::IntPower(
                                  static_cast<std::size_t>(2), dimension - 1)>,
      2 * dimension>
  FaceVertexIndices() {
    using ReturnType =
        std::array<std::array<std::size_t,
                              beziermanipulation::utils::algorithms::IntPower(
                                  static_cast<std::size_t>(2), dimension - 1)>,
                   2 * dimension>;
    if constexpr (dimension == 2) {
      return ReturnType{0, 1, 1, 2, 3, 2, 0, 3};
    } else if constexpr (dimension == 3) {
      // Must be 3D (@todo implement)
      static_assert(dimension == 2, "Not Implemented");
      return ReturnType{};
    }
  }
};
}  // namespace beziermanipulation::utils::algorithms

#endif  // UTILS_TYPE_TRAITS_HYPERCUBE_HPP
