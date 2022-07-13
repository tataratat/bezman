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
#include <cmath>
#include <vector>

#include "bezman/src/rational_bezier_spline.hpp"

using namespace bezman;

namespace bezman::tests::rational_splines_test {

class RationalSplineTestSuite : public ::testing::Test {
 public:
  using Point2D = Point<2, double>;

 protected:
  void SetUp() override {}

  // Provide some data to form splines.
  std::vector<Point2D> ctps{Point2D{1., 0.}, Point2D{1., 1.}, Point2D{0., 1.}};
  std::vector<double> weights{1., 1 / std::sqrt(2), 1.};
  std::vector<Point2D> ctps2{Point2D{0., 0.}, Point2D{0., 1.}, Point2D{1., 1.}};
  std::vector<double> weights2{1., 1., 1.};
  std::array<std::size_t, 1> degree{2};
  std::array<std::size_t, 1> new_degree{3};

  // Testing Point
  Point2D testing_point{2., 2.};
  double testing_weight{3.};
};

// Test Constructors and equal operations
TEST_F(RationalSplineTestSuite, Constructors) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  EXPECT_NO_FATAL_FAILURE(RationalBezier());
  EXPECT_NO_FATAL_FAILURE(RationalBezier(degree));
  EXPECT_NO_FATAL_FAILURE(RationalBezier(degree, ctps, weights));
}

// Getters and Setters
TEST_F(RationalSplineTestSuite, GettersSetters) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  RationalBezier circular_arc(degree, ctps, weights);
  const RationalBezier kcirular_arc(degree, ctps, weights);
  std::array<std::size_t, 1> index_a{0};

  // Get Degrees
  EXPECT_EQ(circular_arc.GetDegrees(), degree);

  // Access Control Points (+ weighted ones) and Weights
  // Weighted CTPS
  EXPECT_EQ(circular_arc.WeightedControlPoint(1), ctps[1] * weights[1]);
  EXPECT_EQ(kcirular_arc.WeightedControlPoint(1), ctps[1] * weights[1]);
  circular_arc.WeightedControlPoint(1) = testing_point;
  EXPECT_EQ(circular_arc.WeightedControlPoint(1), testing_point);
  EXPECT_EQ(circular_arc.WeightedControlPoint(index_a), ctps[0] * weights[0]);
  EXPECT_EQ(kcirular_arc.WeightedControlPoint(index_a), ctps[0] * weights[0]);
  circular_arc.WeightedControlPoint(index_a) = testing_point;
  EXPECT_EQ(circular_arc.WeightedControlPoint(index_a), testing_point);

  circular_arc.WeightedControlPoint(1) = ctps[1] * weights[1];  // Reset arc
  circular_arc.WeightedControlPoint(index_a) = ctps[0] * weights[0];

  // Weights
  EXPECT_EQ(circular_arc.Weight(1), weights[1]);
  EXPECT_EQ(kcirular_arc.Weight(1), weights[1]);
  circular_arc.Weight(1) = testing_weight;
  EXPECT_EQ(circular_arc.Weight(1), testing_weight);
  EXPECT_EQ(circular_arc.Weight(index_a), weights[0]);
  EXPECT_EQ(kcirular_arc.Weight(index_a), weights[0]);
  circular_arc.Weight(index_a) = testing_weight;
  EXPECT_EQ(circular_arc.Weight(index_a), testing_weight);

  circular_arc.Weight(1) = weights[1];
  circular_arc.Weight(index_a) = weights[0];

  // "clean" CTPS
  EXPECT_EQ(circular_arc.ControlPoint(1), ctps[1]);
  EXPECT_EQ(kcirular_arc.ControlPoint(1), ctps[1]);
  EXPECT_EQ(circular_arc.ControlPoint(index_a), ctps[0]);
  EXPECT_EQ(kcirular_arc.ControlPoint(index_a), ctps[0]);

  // Update the degrees
  EXPECT_NO_FATAL_FAILURE(circular_arc.UpdateDegrees(new_degree));
  EXPECT_EQ(circular_arc.GetDegrees(), new_degree);
}

TEST_F(RationalSplineTestSuite, Evaluation) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  RationalBezier circular_arc(degree, ctps, weights);
  Point2D mid_point{std::cos(std::acos(-1) * .25),
                    std::sin(std::acos(-1) * .25)};
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(.5)[0], mid_point[0]);
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(.5)[1], mid_point[1]);
  // Test if symmetric around the x=y plane
  double t{0.15};
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(t)[1],
                  circular_arc.Evaluate(1. - t)[0]);
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(t)[0],
                  circular_arc.Evaluate(1. - t)[1]);
}

TEST_F(RationalSplineTestSuite, ForwardEvaluation) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  RationalBezier circular_arc(degree, ctps, weights);
  Point2D mid_point{std::cos(std::acos(-1) * .25),
                    std::sin(std::acos(-1) * .25)};
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(.5)[0], mid_point[0]);
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(.5)[1], mid_point[1]);
  // Test if symmetric around the x=y plane
  double t{0.15};
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(t)[1],
                  circular_arc.Evaluate(1. - t)[0]);
  EXPECT_FLOAT_EQ(circular_arc.Evaluate(t)[0],
                  circular_arc.Evaluate(1. - t)[1]);
  EXPECT_FLOAT_EQ(circular_arc.ForwardEvaluate(t)[1],
                  circular_arc.ForwardEvaluate(1. - t)[0]);
  EXPECT_FLOAT_EQ(circular_arc.ForwardEvaluate(t)[0],
                  circular_arc.ForwardEvaluate(1. - t)[1]);
}

}  // namespace bezman::tests::rational_splines_test