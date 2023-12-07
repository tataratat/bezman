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
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;
  using Point1D = Point<1, double>;

 protected:
  void SetUp() override {}

  // Provide some data to form splines.
  std::vector<Point2D> ctps{Point2D{1., 0.}, Point2D{1., 1.}, Point2D{0., 1.}};
  std::vector<double> weights{1., 1 / std::sqrt(2), 1.};
  std::vector<Point2D> ctps2{Point2D{0., 0.}, Point2D{0., 1.}, Point2D{1., 1.}};
  std::vector<double> weights2{1., 1., 1.};
  std::array<std::size_t, 1> degree{2};
  std::array<std::size_t, 1> new_degree{3};

  // Higher Dimensional Examples
  std::vector<Point2D> ctps_2D{
      Point2D{1., 0.}, Point2D{1., 1.}, Point2D{0., 1.},
      Point2D{2., 0.}, Point2D{2., 2.}, Point2D{0., 2.},
  };
  std::vector<double> weights_2D{1., 1 / std::sqrt(2), 1.,
                                 1., 1 / std::sqrt(2), 1.};
  std::vector<Point2D> ctps2_2D{
      Point2D{0., 0.}, Point2D{1., 0.}, Point2D{2., 0.},
      Point2D{0., 1.}, Point2D{1., 1.}, Point2D{2., 1.},
      Point2D{0., 2.}, Point2D{1., 2.}, Point2D{2., 2.}};
  std::vector<double> weights2_2D{1., 1., 1., 1., 1., 1., 1., 1., 1.};
  std::array<std::size_t, 2> degrees_2D{2, 1};
  std::array<std::size_t, 2> degrees2_2D{2, 2};

  // Testing Point
  Point2D testing_point{2., 2.};
  double testing_weight{3.};

  // Testing a derivative with a scalar bezier
  std::vector<double> ctps_deriv_test{3., 4., 5.};
  std::vector<double> weights_deriv_test{1., .5, 1.};

  double AnalyticalSolutionToRspline(const double x) {
    return (1. + 2. * x - 2 * x * x) / std::pow(1 - x + x * x, 2);
  }

  // Create a randomized spline
  template <std::size_t para_dim>
  auto CreateRandomSpline(const std::array<std::size_t, para_dim> &degrees) {
    RationalBezierSpline<para_dim, Point3D, double> randomSpline{degrees};
    for (std::size_t i_ctps{}; i_ctps < randomSpline.GetNumberOfControlPoints();
         i_ctps++) {
      for (std::size_t i_dim{}; i_dim < 3; i_dim++) {
        randomSpline.GetWeightedControlPoints()[i_ctps][i_dim] =
            static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
      }
      randomSpline.GetWeights()[i_ctps] =
          static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }
    return randomSpline;
  }
};

// Test Constructors and equal operations
TEST_F(RationalSplineTestSuite, Constructors) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  EXPECT_NO_FATAL_FAILURE(RationalBezier());
  EXPECT_NO_FATAL_FAILURE(RationalBezier{degree});
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

// Test Basis Function evaluation
TEST_F(RationalSplineTestSuite, TestBasisFunctions) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  RationalBezier circular_arc(degree, ctps, weights);
  // Create random evaluation point
  const double xx{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  const Point1D x{xx};
  // Retrieve basis functions (check overloads)
  EXPECT_EQ(circular_arc.BasisFunctions(xx), circular_arc.BasisFunctions(x));
  EXPECT_EQ(circular_arc.BasisFunctionContributions(xx),
            circular_arc.BasisFunctionContributions(x));

  // Compare to analytical solution
  double weighted_basis_function_sum{};

  // Test Polynomial Basis Functions
  const auto basis_functions = circular_arc.BasisFunctionContributions(x);
  for (std::size_t i_basis{}; i_basis < degree[0] + 1; i_basis++) {
    const double analytical_solution_bernstein_pol =
        utils::FastBinomialCoefficient::choose(degree[0], i_basis) *
        std::pow(x[0], i_basis) * std::pow(1. - x[0], degree[0] - i_basis);
    weighted_basis_function_sum +=
        analytical_solution_bernstein_pol * circular_arc.Weight(i_basis);
    EXPECT_FLOAT_EQ(basis_functions[0][i_basis],
                    analytical_solution_bernstein_pol);
  }

  // Test Weighted Basis Functions
  const auto weighted_basis_function = circular_arc.BasisFunctions(xx);
  const auto non_weighted_basis_function =
      circular_arc.UnweightedBasisFunctions(xx);

  for (std::size_t i_basis{}; i_basis < degree[0] + 1; i_basis++) {
    const double analytical_solution_bernstein_pol =
        utils::FastBinomialCoefficient::choose(degree[0], i_basis) *
        std::pow(x[0], i_basis) * std::pow(1. - x[0], degree[0] - i_basis);
    EXPECT_DOUBLE_EQ(weighted_basis_function[i_basis],
                     analytical_solution_bernstein_pol *
                         circular_arc.GetWeights()[i_basis] /
                         weighted_basis_function_sum);
    EXPECT_DOUBLE_EQ(
        non_weighted_basis_function[i_basis],
        analytical_solution_bernstein_pol / weighted_basis_function_sum);
  }
}

// Test High Dimensional Evaluation
TEST_F(RationalSplineTestSuite, HighDimEvaluation) {
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  RationalBezierSurface circle_segment{degrees_2D, ctps_2D, weights_2D};
  // Evaluate
  Point2D mid_point{std::cos(std::acos(-1) * .25) * 1.5,
                    std::sin(std::acos(-1) * .25) * 1.5};
  EXPECT_FLOAT_EQ(circle_segment.Evaluate(.5, .5)[0], mid_point[0]);
  EXPECT_FLOAT_EQ(circle_segment.Evaluate(.5, .5)[1], mid_point[1]);

  // Test Basis Functions
  RationalBezierSurface pseudo_polynomial_rectangle{degrees2_2D, ctps2_2D,
                                                    weights2_2D};
  const auto degrees = pseudo_polynomial_rectangle.GetDegrees();
  // Create random evaluation point
  const double xx{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  const double xy{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
  const Point2D x{xx, xy};
  // Retrieve basis functions
  const auto basis_functions =
      pseudo_polynomial_rectangle.BasisFunctions(xx, xy);
  const auto pol_basis_functions =
      pseudo_polynomial_rectangle.BasisFunctions(xx, xy);

  // Compare to analytical solution
  for (std::size_t pdim{}; pdim < 2; pdim++) {
    for (std::size_t i_basis{}; i_basis < degrees[pdim] + 1; i_basis++) {
      EXPECT_FLOAT_EQ(basis_functions[i_basis], pol_basis_functions[i_basis]);
    }
  }
}

// Test Refinement
TEST_F(RationalSplineTestSuite, RefinementTest) {
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  RationalBezierSurface circle_segment{degrees_2D, ctps_2D, weights_2D};
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  RationalBezier circular_arc(degree, ctps, weights);
  const RationalBezierSurface kcircle_segment{degrees_2D, ctps_2D, weights_2D};
  const RationalBezier kcircular_arc(degree, ctps, weights);

  // Test specification
  const std::size_t n_refinement_level{3};
  const std::size_t n_sample_points{10};

  // Refinement
  for (std::size_t i_refinement{0}; i_refinement < n_refinement_level;
       i_refinement++) {
    circle_segment.OrderElevateAlongParametricDimension(0);
    circle_segment.OrderElevateAlongParametricDimension(1);
    circular_arc.OrderElevateAlongParametricDimension(0);
  }

  // Sampling
  for (std::size_t i_sample_x{0}; i_sample_x < n_sample_points; i_sample_x++) {
    const double xx{static_cast<double>(rand()) /
                    static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(circular_arc.Evaluate(xx)[0],
                    kcircular_arc.Evaluate(xx)[0]);
    for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points;
         i_sample_y++) {
      const double xy{static_cast<double>(rand()) /
                      static_cast<double>(RAND_MAX)};
      EXPECT_FLOAT_EQ(circle_segment.Evaluate(xx, xy)[0],
                      kcircle_segment.Evaluate(xx, xy)[0]);
      EXPECT_FLOAT_EQ(circle_segment.Evaluate(xx, xy)[1],
                      kcircle_segment.Evaluate(xx, xy)[1]);
    }
  }
}

// Test Derivatives of rational Beziers
TEST_F(RationalSplineTestSuite, DerivativeTesting) {
  using RationalBezier = RationalBezierSpline<1, double, double>;
  RationalBezier spline{degree, ctps_deriv_test, weights_deriv_test};
  const auto derivative =
      spline.DerivativeWRTParametricDimension(std::array<std::size_t, 1>{1});
  // Sample
  const std::size_t n_sample_points{10};
  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(derivative.Evaluate(x), AnalyticalSolutionToRspline(x));
  }
}

// Test Derivatives of rational Beziers
TEST_F(RationalSplineTestSuite, DerivativeEvaluation) {
  constexpr std::size_t n_tests{3}, para_dims{2};
  const std::size_t n_test_evaluations{5};
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    const auto degrees = std::array<std::size_t, para_dims>{3, 3};
    const auto derivs = std::array<std::array<std::size_t, para_dims>, n_tests>{
        {{2, 0}, {1, 1}, {1, 2}}}[i_test];
    const auto random_spline = CreateRandomSpline(degrees);
    auto random_spline_deriv = random_spline;

    random_spline_deriv =
        random_spline_deriv.DerivativeWRTParametricDimension(derivs);

    // Evaluate for comparison
    for (std::size_t j_test{}; j_test < n_test_evaluations; j_test++) {
      // Create evaluation points
      const double x{static_cast<double>(rand()) /
                     static_cast<double>(RAND_MAX)};
      const double y{static_cast<double>(rand()) /
                     static_cast<double>(RAND_MAX)};
      Point2D eval_point{x, y};
      const auto result = random_spline.EvaluateDerivative(
          // Evaluation Point
          eval_point,
          // derivatives
          derivs);
      const auto result_2 = random_spline_deriv.Evaluate(eval_point);
      // Compare results dimensions
      for (std::size_t i_dim{}; i_dim < 3; i_dim++) {
        EXPECT_FLOAT_EQ(result[i_dim], result_2[i_dim]);
      }
    }
  }
}

// Test Derivatives of rational Beziers
TEST_F(RationalSplineTestSuite, BasisFunctionDerivativeEvaluation) {
  constexpr std::size_t n_tests{4}, para_dims{2};
  const std::size_t n_test_evaluations{5};
  for (std::size_t i_test{}; i_test < n_tests; i_test++) {
    const auto degrees = std::array<std::size_t, para_dims>{3, 3};
    const auto derivs = std::array<std::array<std::size_t, para_dims>, n_tests>{
        {{2, 0}, {1, 1}, {1, 2}, {1, 3}}}[i_test];
    const auto random_spline = CreateRandomSpline(degrees);

    // Evaluate for comparison
    for (std::size_t j_test{}; j_test < n_test_evaluations; j_test++) {
      // Create evaluation points
      const double x{static_cast<double>(rand()) /
                     static_cast<double>(RAND_MAX)};
      const double y{static_cast<double>(rand()) /
                     static_cast<double>(RAND_MAX)};
      Point2D eval_point{x, y};
      const auto result = random_spline.EvaluateDerivative(
          // Evaluation Point
          eval_point,
          // derivatives
          derivs);
      const auto basis_function_derivatives =
          random_spline.BasisFunctionsDerivatives(
              // Evaluation Point
              eval_point,
              // derivatives
              derivs);
      Point3D sum{};
      for (std::size_t i_basis{};
           i_basis < random_spline.GetNumberOfControlPoints(); i_basis++) {
        sum += basis_function_derivatives[i_basis] *
               random_spline.GetWeightedControlPoints()[i_basis] /
               random_spline.GetWeights()[i_basis];
      };
      // Compare results dimensions
      for (std::size_t i_dim{}; i_dim < 3; i_dim++) {
        EXPECT_FLOAT_EQ(result[i_dim], sum[i_dim]);
      }
    }
  }
}

// Test Multiplication
TEST_F(RationalSplineTestSuite, MultiplicationTest) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  RationalBezier circular_arc_v(degree, ctps, weights);
  // Create random evaluation point
  const std::size_t n_sample_points{10};

  const double scalar_v{static_cast<double>(rand()) /
                        static_cast<double>(RAND_MAX)};
  circular_arc_v *= scalar_v;
  const auto circular_multiplied = circular_arc * scalar_v;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(circular_multiplied.Evaluate(x)[0],
                    circular_arc.Evaluate(x)[0] * scalar_v);
    EXPECT_FLOAT_EQ(circular_multiplied.Evaluate(x)[1],
                    circular_arc.Evaluate(x)[1] * scalar_v);
    EXPECT_FLOAT_EQ(circular_multiplied.Evaluate(x)[0],
                    circular_arc_v.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(circular_multiplied.Evaluate(x)[1],
                    circular_arc_v.Evaluate(x)[1]);
  }
}

// Test Addition
TEST_F(RationalSplineTestSuite, AdditionTest) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  RationalBezier circular_arc_v(degree, ctps, weights);
  // Create random evaluation point
  const std::size_t n_sample_points{10};

  const Point2D shifting_vector{
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX),
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};

  circular_arc_v += shifting_vector;
  const auto circular_modified = circular_arc + shifting_vector;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc.Evaluate(x)[0] + shifting_vector[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc.Evaluate(x)[1] + shifting_vector[1]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc_v.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc_v.Evaluate(x)[1]);
  }
}

// Test Substraction
TEST_F(RationalSplineTestSuite, Substraction) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  RationalBezier circular_arc_v(degree, ctps, weights);
  // Create random evaluation point
  const std::size_t n_sample_points{10};

  const Point2D shifting_vector{
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX),
      static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};

  circular_arc_v -= shifting_vector;
  const auto circular_modified = circular_arc - shifting_vector;
  const auto inverted_spline = -circular_arc;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc.Evaluate(x)[0] - shifting_vector[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc.Evaluate(x)[1] - shifting_vector[1]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc_v.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc_v.Evaluate(x)[1]);
    EXPECT_FLOAT_EQ(inverted_spline.Evaluate(x)[0],
                    -circular_arc.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(inverted_spline.Evaluate(x)[1],
                    -circular_arc.Evaluate(x)[1]);
  }
}

// Test Addition
TEST_F(RationalSplineTestSuite, SplineAdditionTest) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  RationalBezier circular_arc_v(degree, ctps, weights);
  const RationalBezier secondary_spline{degree, ctps2, weights2};

  // Create random evaluation point
  const std::size_t n_sample_points{10};

  // Adding Splines with the same weights does not elevate degree
  EXPECT_EQ((circular_arc + circular_arc).GetDegrees(),
            circular_arc.GetDegrees());

  circular_arc_v += secondary_spline;
  const auto circular_modified = circular_arc + secondary_spline;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(
        circular_modified.Evaluate(x)[0],
        circular_arc.Evaluate(x)[0] + secondary_spline.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(
        circular_modified.Evaluate(x)[1],
        circular_arc.Evaluate(x)[1] + secondary_spline.Evaluate(x)[1]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc_v.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc_v.Evaluate(x)[1]);
  }
}

// Test Substraction of Splines
TEST_F(RationalSplineTestSuite, SplineSubstractionTest) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  RationalBezier circular_arc_v(degree, ctps, weights);
  const RationalBezier secondary_spline{degree, ctps2, weights2};

  // Create random evaluation point
  const std::size_t n_sample_points{10};

  // Adding Splines with the same weights does not elevate degree
  EXPECT_EQ((circular_arc - circular_arc).GetDegrees(),
            circular_arc.GetDegrees());

  circular_arc_v -= secondary_spline;
  const auto circular_modified = circular_arc - secondary_spline;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(
        circular_modified.Evaluate(x)[0],
        circular_arc.Evaluate(x)[0] - secondary_spline.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(
        circular_modified.Evaluate(x)[1],
        circular_arc.Evaluate(x)[1] - secondary_spline.Evaluate(x)[1]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[0],
                    circular_arc_v.Evaluate(x)[0]);
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x)[1],
                    circular_arc_v.Evaluate(x)[1]);
  }
}

// Test Multiplication of Splines
TEST_F(RationalSplineTestSuite, SplineMultiplicationTest) {
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  const RationalBezier secondary_spline{degree, ctps2, weights2};

  // Create random evaluation point
  const std::size_t n_sample_points{10};

  const auto circular_modified = circular_arc * secondary_spline;

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    EXPECT_FLOAT_EQ(circular_modified.Evaluate(x),
                    circular_arc.Evaluate(x) * secondary_spline.Evaluate(x));
  }
}

// Test Unit Cube Checks
TEST_F(RationalSplineTestSuite, CubeCheckSplineTest) {
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  const RationalBezierSurface circle_segment{degrees_2D, ctps_2D, weights_2D};
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  EXPECT_FLOAT_EQ(circle_segment.MaximumCorner()[0], Point2D(2., 2.)[0]);
  EXPECT_FLOAT_EQ(circle_segment.MaximumCorner()[1], Point2D(2., 2.)[1]);
  EXPECT_FLOAT_EQ(circle_segment.MinimumCorner()[0], Point2D(0., 0.)[0]);
  EXPECT_FLOAT_EQ(circle_segment.MinimumCorner()[1], Point2D(0., 0.)[1]);
  EXPECT_FALSE(circle_segment.FitsIntoUnitCube());
  EXPECT_TRUE(circular_arc.FitsIntoUnitCube());
}

// Test Composition of (Rational) splines
TEST_F(RationalSplineTestSuite, SplineComposition) {
  using RationalBezierSurface = RationalBezierSpline<2, Point2D, double>;
  RationalBezierSurface circle_segment_large =
      RationalBezierSurface(degrees_2D, ctps_2D, weights_2D);
  RationalBezierSurface circle_segment_small = circle_segment_large * 0.5;
  using RationalBezier = RationalBezierSpline<1, Point2D, double>;
  using BezierSurface = BezierSpline<2, Point2D, double>;
  const RationalBezier circular_arc(degree, ctps, weights);
  const BezierSurface rectangle{
      // degrees
      std::array<std::size_t, 2>{1, 1},
      // CTPS
      std::vector<Point2D>{Point2D{0.5, 0.}, Point2D{1., 0.3}, Point2D{0., 0.3},
                           Point2D{0.5, 1.}}};

  // Create random evaluation point
  const std::size_t n_sample_points{5};

  /// Perform Compositions
  // Polynomial (Rational) [1D -> 2D]
  const auto circular_composed = rectangle.Compose(circular_arc);
  // Polynomial (Rational) [2D -> 2D]
  const auto circle_segment_in_rectangle =
      rectangle.Compose(circle_segment_small);
  // Rational (Rational) [2D -> 2D]
  const auto circle_circle_composition =
      circle_segment_large.Compose(circle_segment_small);
  // Rational (Polynomial)
  const auto rectangle_in_circle_composition =
      circle_segment_large.Compose(rectangle);

  for (std::size_t i_sample_y{0}; i_sample_y < n_sample_points; i_sample_y++) {
    // Sample points
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const double y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    // Polynomial (Rational) [1D -> 2D]
    EXPECT_FLOAT_EQ(circular_composed.Evaluate(x)[0],
                    rectangle.Evaluate(circular_arc.Evaluate(x))[0]);
    EXPECT_FLOAT_EQ(circular_composed.Evaluate(x)[1],
                    rectangle.Evaluate(circular_arc.Evaluate(x))[1]);
    // Polynomial (Rational) [2D -> 2D]
    EXPECT_FLOAT_EQ(circle_segment_in_rectangle.Evaluate(x, y)[0],
                    rectangle.Evaluate(circle_segment_small.Evaluate(x, y))[0]);
    EXPECT_FLOAT_EQ(circle_segment_in_rectangle.Evaluate(x, y)[1],
                    rectangle.Evaluate(circle_segment_small.Evaluate(x, y))[1]);
    // Rational (Rational) [2D -> 2D]
    EXPECT_FLOAT_EQ(
        circle_circle_composition.Evaluate(x, y)[0],
        circle_segment_large.Evaluate(circle_segment_small.Evaluate(x, y))[0]);
    EXPECT_FLOAT_EQ(
        circle_circle_composition.Evaluate(x, y)[1],
        circle_segment_large.Evaluate(circle_segment_small.Evaluate(x, y))[1]);
    // Rational (Rational) [2D -> 2D]
    EXPECT_FLOAT_EQ(rectangle_in_circle_composition.Evaluate(x, y)[0],
                    circle_segment_large.Evaluate(rectangle.Evaluate(x, y))[0]);
    EXPECT_FLOAT_EQ(rectangle_in_circle_composition.Evaluate(x, y)[1],
                    circle_segment_large.Evaluate(rectangle.Evaluate(x, y))[1]);
  }
}

}  // namespace bezman::tests::rational_splines_test
