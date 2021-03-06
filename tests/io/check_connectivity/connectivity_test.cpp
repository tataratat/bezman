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
   *               0 --(0)-- 1
   */
  using Point2D = ::bezman::Point<2ul, double>;

 public:
  // Face center points
  std::vector<Point2D> face_center_points{
      // First Patch
      Point2D{0.5, 1.}, Point2D{1., 1.5}, Point2D{0.5, 2.}, Point2D{0., 1.5},
      // Second Patch
      Point2D{1.5, 1.}, Point2D{2., 1.5}, Point2D{1.5, 2.}, Point2D{1., 1.5},
      // Third Patch
      Point2D{1.5, 0.}, Point2D{2., 0.5}, Point2D{1.5, 1.}, Point2D{1., 0.5}};

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

  // List of opposite faces
  std::array<std::size_t, 4> opposite_faces{2, 3, 0, 1};

  // Metric
  Point2D metric{1.1, 1.};

  // Unique List
  std::vector<std::size_t> uniquelist{1, 5, 4, 0, 6, 9, 8, 5, 3, 7, 6, 2};
};

/*
 * Check if the connectivity is calculated correctly
 */
TEST_F(ConnectivityCheckSuite, ConnectivityTest) {
  // Calculate the connectivity
  const auto connectivity =
      algorithms::FindConnectivity(face_center_points, metric, opposite_faces);

  EXPECT_EQ(connectivity, expected_connectivity);
}

/*
 * Check if the uniquifying methods for single points works as well
 */
TEST_F(ConnectivityCheckSuite, UniquifyPointsTest) {
  // Calculate the connectivity
  const auto connectivity =
      algorithms::IndexUniquePointList(face_center_points, metric);

  EXPECT_EQ(connectivity, uniquelist);
}

}  // namespace bezman::tests::io::connectivity_check