/*
 * This Example illustrates, how the BezierManipulation suite can be used to
 * build parametrized microstructures. It is also used to show the capabilities
 * of deriving the structure with respect to the design variables, using
 * algorithmic differentiation.
 *
 * Implementationwise, the code is structured as follows. There are three
 * classes describing the geometry and the parametrization. The MicortileExample
 * class, describes the parametrized Inner-spline, that is later inserted into
 * the structure. It contains a function, that creates a microtile, based on a
 * set of parameters. As the parameters are based on a function evaluation, it
 * also provides the points, where this funciton is to be evaluated (in the case
 * of a crosstile, these are the tips of the crosstile-branches). The
 * DeformationFunctionExample class provides the outer geometry. As it is formed
 * through a group of Bezier splines, it also holds segmentation information,
 * i.e., how many splines are in each parametric dimension. This could also be
 * generalized with a function, returning local coordinates based on the current
 * index. The ValueFunction class finally provides the parametrization function,
 * that returns the local Parameters, based on the position of the spline within
 * the structure and a set of superordinate parameters. These three functions
 * are then inserted as Template Arguments in the Microstructure-Generator
 * class, which performes the compositions and derivations.
 *
 * The results are exported in xml format and can be plotted using the python
 * scripts with gustav
 *
 */

#include <algorithm>
#include <cmath>
#include <vector>

#include "bezierManipulation/src/bezier_spline_group.hpp"
#include "bezierManipulation/src/utils/computational_derivation/algo_diff_type.hpp"
#include "bezierManipulation/src/utils/export.hpp"
#include "cross_tile_3d.hpp"
#include "microstructure_generator.hpp"
#include "ring_segments_3d.hpp"

using namespace beziermanipulation;

/**
 * @brief Value Field Function defined in the parametric domain of the
 * deformation function
 *
 * This function is based on Algorithmic differentiation types, which
 calculate
 * the derivative of the value field function at the same time.
 *
 */
class ValueFieldExample {
 private:
  using ADT = utils::computational_derivation::AlgoDiffType<double>;
  using Point3D = Point<3, double>;

 public:
  // Number of super Parameters
  static constexpr const std::size_t kNumberOfSuperParameters{12};

  // Some super parameters
  const std::array<ADT, kNumberOfSuperParameters> kSuperControlPoints{
      ADT(0.1, kNumberOfSuperParameters, 0),
      ADT(0.2, kNumberOfSuperParameters, 1),
      ADT(0.1, kNumberOfSuperParameters, 2),
      ADT(0.1, kNumberOfSuperParameters, 3),
      ADT(0.1, kNumberOfSuperParameters, 4),
      ADT(0.3, kNumberOfSuperParameters, 5),
      ADT(0.1, kNumberOfSuperParameters, 6),
      ADT(0.2, kNumberOfSuperParameters, 7),
      ADT(0.1, kNumberOfSuperParameters, 8),
      ADT(0.1, kNumberOfSuperParameters, 9),
      ADT(0.1, kNumberOfSuperParameters, 10),
      ADT(0.3, kNumberOfSuperParameters, 11)};

  /**
   * @brief Evaluate the Value field function
   *
   * Here, the value field function are the basis functions of a B-Spline
   * interpolating between the SuperControlPoints. The B-Spline has a knot in
   * the center, meaning that the basis functions have only local support.
   */
  template <std::size_t array_size>
  std::array<ADT, array_size> Evaluate(
      const std::array<Point3D, array_size>& evaluation_points) const {
    // Initialize return value
    std::array<ADT, array_size> evaluations{};
    // Assign based on simple linear B-Spline basis
    for (std::size_t i_point{}; i_point < array_size; i_point++) {
      const double &x{evaluation_points[i_point][0]},
          &y{evaluation_points[i_point][1]}, &z{evaluation_points[i_point][2]};
      evaluations[i_point] =
          (x < 0.5 ? (1 - 2. * x) : 0.) * (1 - y) * (1 - z) *
              kSuperControlPoints[0] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * (1 - y) * (1 - z) *
              kSuperControlPoints[1] +
          (x < 0.5 ? 0. : 2 * x - 1.) * (1 - y) * (1 - z) *
              kSuperControlPoints[2] +
          (x < 0.5 ? (1 - 2. * x) : 0.) * y * (1 - z) * kSuperControlPoints[3] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * y * (1 - z) *
              kSuperControlPoints[4] +
          (x < 0.5 ? 0. : 2 * x - 1) * y * (1 - z) * kSuperControlPoints[5] +
          (x < 0.5 ? (1 - 2. * x) : 0.) * (1 - y) * z * kSuperControlPoints[6] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * (1 - y) * z *
              kSuperControlPoints[7] +
          (x < 0.5 ? 0. : 2 * x - 1.) * (1 - y) * z * kSuperControlPoints[8] +
          (x < 0.5 ? (1 - 2. * x) : 0.) * y * z * kSuperControlPoints[9] +
          (x < 0.5 ? (2. * x) : (2. - 2. * x)) * y * z *
              kSuperControlPoints[10] +
          (x < 0.5 ? 0. : 2 * x - 1) * y * z * kSuperControlPoints[11];
    }
    return evaluations;
  }
};

int main() {
#ifdef ENABLE_OPEN_MP_PARALLEL
  omp_set_num_threads(4);
#endif
  // Initialize generator
  auto micro_structure_generator =
      MicrostructureGenerator<CrossTile3D, RingSegments3D, ValueFieldExample>{};
  // Modify the deformation function
  micro_structure_generator.deformation_function_generator.SetNumberOfSegments(
      std::array<int, 3>{15, 15, 15});

  // Construct the composition
  const auto test_composition =
      micro_structure_generator.ComposeMicrostructureAndDerivatives();
  std::cout << "Starting the file export" << std::endl;
  // Export the MS and its derivatives
  utils::Export::GuessByExtension(test_composition[0],
                                  "composed_microstructure.itd");

#ifdef ENABLE_OPEN_MP_PARALLEL
#pragma omp parallel for
#endif
  for (std::size_t i_deriv = std::size_t{};
       i_deriv < ValueFieldExample::kNumberOfSuperParameters; i_deriv++) {
    utils::Export::GuessByExtension(test_composition[i_deriv + 1],
                                    std::string("composed_microstructure_") +
                                        std::to_string(i_deriv) +
                                        std::string(".itd"));
  }

  return 0;
}
