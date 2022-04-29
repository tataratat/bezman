#include <gtest/gtest.h>

#include <array>
#include <iostream>
#include <vector>

#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/uniquify/point_uniquifier.hpp"

using namespace beziermanipulation::utils;

namespace beziermanipulation::tests::io::connectivity_check {

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
  using Point2D = ::beziermanipulation::Point<2ul, double>;

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
  std::vector<std::array<int, 4>> expected_connectivity{
      {-1, 1, -1, -1}, {2, -1, -1, 0}, {-1, -1, 1, -1}};

  // List of opposite faces
  std::array<std::size_t, 4> opposite_faces{2, 3, 0, 1};

  // Metric
  Point2D metric{0.25, 1.};
};

/*
 * Check if the connectivity is calculated correctly
 */
TEST_F(ConnectivityCheckSuite, ConnectivityTest) {
  // Calculate the connectivity
  const auto connectivity =
      uniquify::FindConnectivity(face_center_points, metric, opposite_faces);

  EXPECT_EQ(connectivity, expected_connectivity);
}

}  // namespace beziermanipulation::tests::io::connectivity_check