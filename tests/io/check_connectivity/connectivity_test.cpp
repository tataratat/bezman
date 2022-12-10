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

#include <gtest/gtest.h>

#include <array>
#include <iostream>
#include <vector>

#include "bezman/src/point.hpp"
#include "bezman/src/utils/algorithms/point_uniquifier.hpp"

using namespace bezman::utils;

namespace bezman::tests::io::connectivity_check {

class ConnectivityCheckSuite : public ::testing::Test {
  /*
   * Arrangement:
   *
   *  3 --(2)-- 2  3 --(2)-- 2
   *  |         |  |         |
   * (3)   0   (1)(3)   1   (1)
   *  |         |  |         |
   *  0 --(0)-- 1  0 --(0)-- 1
   *               3 --(2)-- 2
   *               |         |
   *              (3)   2   (1)
   *               |         |
   * origin        0 --(0)-- 1
   */
  using Point2D = ::bezman::Point<2ul, double>;
  using Point3D = ::bezman::Point<3ul, double>;

 public:
  // Face center points
  std::vector<Point2D> corner_vertices{
      // First Patch
      Point2D{0., 1.}, Point2D{1., 1.}, Point2D{1., 2.}, Point2D{0., 2.},
      // Second Patch
      Point2D{1., 1.}, Point2D{2., 1.}, Point2D{2., 2.}, Point2D{1., 2.},
      // Third Patch
      Point2D{1., 0.}, Point2D{2., 0.}, Point2D{2., 1.}, Point2D{1., 1.}};

  // Face center points embdedded
  std::vector<Point3D> corner_vertices_3d{
      // First Patch
      Point3D{0., 1., 2.}, Point3D{1., 1., 2.}, Point3D{1., 2., 2.},
      Point3D{0., 2., 2.},
      // Second Patch
      Point3D{1., 1., 2.}, Point3D{2., 1., 2.}, Point3D{2., 2., 2.},
      Point3D{1., 2., 2.},
      // Third Patch
      Point3D{1., 0., 2.}, Point3D{2., 0., 2.}, Point3D{2., 1., 2.},
      Point3D{1., 1., 2.}};

  // Expected Connectivity
  std::vector<std::array<std::size_t, 4>> expected_connectivity{
      // First element
      {static_cast<std::size_t>(-1),   // 0
       1,                              // 1
       static_cast<std::size_t>(-1),   // 2
       static_cast<std::size_t>(-1)},  // 3
      // Second element
      {2,                             // 0
       static_cast<std::size_t>(-1),  // 1
       static_cast<std::size_t>(-1),  // 2
       0},                            // 3
      // Third element
      {static_cast<std::size_t>(-1),  // 0
       static_cast<std::size_t>(-1),  // 1
       1,                             // 2
       static_cast<std::size_t>(-1)}  // 3
  };

  // Metric 2D
  Point2D metric{1.1, 1.};
  // Metric 3D
  Point3D metric3d{1.1, 1., 0.1};

  // Unique List
  std::vector<std::size_t> uniquelist{0, 3, 5, 2, 3, 6, 7, 5, 1, 4, 6, 3};
};

/*
 * Check if the connectivity is calculated correctly
 */
TEST_F(ConnectivityCheckSuite, ConnectivityTest) {
  // Calculate the connectivity
  const auto connectivity =
      algorithms::FindConnectivityFromCorners<2>(corner_vertices, metric, 1e-5);

  EXPECT_EQ(connectivity, expected_connectivity);
}

/*
 * Check the connectivity for an embedded problem
 */
TEST_F(ConnectivityCheckSuite, ConnectivityTestEmbedded) {
  // Calculate the connectivity
  const auto connectivity = algorithms::FindConnectivityFromCorners<2>(
      corner_vertices_3d, metric3d, 1e-5);

  EXPECT_EQ(connectivity, expected_connectivity);
}

/*
 * Check if the uniquifying methods for single points works as well
 */
TEST_F(ConnectivityCheckSuite, UniquifyPointsTest) {
  // Calculate the connectivity
  const auto connectivity =
      algorithms::IndexUniquePointList(corner_vertices, metric);

  EXPECT_EQ(connectivity, uniquelist);
}

}  // namespace bezman::tests::io::connectivity_check
