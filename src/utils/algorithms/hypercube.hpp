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
    template<std::size_t dimension>
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
        static constexpr std::array<std::size_t, dimension * 2>
        GetOppositeFaces() {
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
        static std::array<std::size_t, beziermanipulation::utils::algorithms::IntPower(
                static_cast<std::size_t>(2), dimension)> VertexIdForDegrees(
                const std::array<std::size_t, dimension> &degrees) {

            using ReturnType = std::array<std::size_t, beziermanipulation::utils::algorithms::IntPower(
                    static_cast<std::size_t>(2), dimension)>;
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
                return ReturnType{
                        0,
                        degrees[0],
                        (degrees[0] + 1) * (degrees[1] + 1) - 1,
                        (degrees[0] + 1) * degrees[1]};
            } else {
                // Must be 3D (@todo implement)
                static_assert(dimension == 2, "Not Implemented");

                /**
                 *
                 *   Here could be a cube with the corners and edges
                 *
                 */

                return std::array<std::array<std::size_t, 4>, dimension * 2>{
                        // bottom: 0 - 1 - 2 - 3                                // vertex
                        0,                                                      //0
                        degrees[0],                                             //1
                        (degrees[0] + 1) * (degrees[1] + 1) -
                        1,                //2
                        (degrees[0] + 1) *
                        degrees[1],                          //3

                        // top: 4 - 5 - 6 - 7
                        degrees[2],                                             //4
                        degrees[2] +
                        degrees[0],                                //5
                        degrees[2] + (degrees[0] + 1) * (degrees[1] + 1) -
                        1,    //6
                        degrees[2] +
                        (degrees[0] + 1) * degrees[1],             //7

                };
            }
        }
            /**
             * @brief Vertex indices (of corners) on faces
             *
             * @return constexpr std::array<std::size_t, 2^(dim-1)> vertices on faces
             */
            static constexpr std::array<
                    std::array<
                            std::size_t,
                            static_cast<std::size_t>(2)
                    >,
                    dimension == 2 ? static_cast<std::size_t>(4)
                                   : static_cast<std::size_t>(12)
            > EdgeVertexIndices() {

                using ReturnType =
                        std::array<std::array<std::size_t,
                                beziermanipulation::utils::algorithms::IntPower(
                                        static_cast<std::size_t>(2),
                                        dimension - 1)>,
                                2 * dimension>;

                static_assert(dimension == 2 || dimension == 3,
                              "Not Implemented");
                if constexpr (dimension == 2) {
                    return ReturnType{0, 1, 1, 2, 3, 2, 0, 3};
                } else if constexpr (dimension == 3) {
                    return ReturnType{
                            // Horizontal edges bottom
                            0, 1, 1, 2, 3, 2, 0, 3,
                            // Horizontal Edges Top
                            4, 5, 5, 6, 7, 6, 4, 7,
                            // Vertical edges
                            0, 4, 1, 5, 2, 6, 3, 7
                    };
                }
            };

    };
}// namespace beziermanipulation::utils::algorithms
#endif  // UTILS_TYPE_TRAITS_HYPERCUBE_HPP
