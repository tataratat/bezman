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

#ifndef EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROSTRUCTURE_GENERATOR_HPP
#define EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROSTRUCTURE_GENERATOR_HPP

#ifdef ENABLE_OPEN_MP_PARALLEL
#include "omp.h"
#endif

/**
 * @brief Generator of the parametrized microstructure
 *
 * Constructs the parametrized microstructure along with its derivatives
 *
 * @tparam Microtile            Microtile Example (Parametrized function)
 * @tparam DeformationFunction  Deformation function as a PolyBezierGroup
 * @tparam ValueField           Value function that is used to determine the
 * local Microtileparameters
 */
template <typename Microtile, typename DeformationFunction, typename ValueField>
class MicrostructureGenerator {
 private:
  // Aliases
  using ADT = utils::computational_differentiation::AlgoDiffType<double>;
  using PointADT3D = Point<3, ADT>;
  using Point3D = Point<3, double>;
  using Bezier = BezierSpline<3, Point3D, double>;
  using PolyBezierGroup = BezierGroup<Bezier>;

  // TransformPoints_ transforms the points provided by
  // Microtile::kEvaluationPoints these pints are transformed by a linear
  // interpolation funciton and set into the Valuefield function to retrieve
  // the local parameters  of the microtile function.
  template <std::size_t number_of_points>
  std::array<Point3D, number_of_points> TransformPoints_(
      const Point3D& cornermin, const Point3D& cornermax,
      const std::array<Point3D, number_of_points>& evaluation_points) const {
    std::array<Point3D, number_of_points> return_value{};
    const auto difference = cornermax - cornermin;
    for (std::size_t i{}; i < number_of_points; i++) {
      return_value[i] =
          cornermin + Point3D{evaluation_points[i][0] * difference[0],
                              evaluation_points[i][1] * difference[1],
                              evaluation_points[i][2] * difference[2]};
    }
    return return_value;
  }

 public:
  // Instance of the value field function, MicrotileGenerator and
  // DeformationFunction, they should be publicly available, to facilitate their
  // modification
  DeformationFunction deformation_function_generator{};
  ValueField value_field{};
  Microtile microtile{};

  // Compose microstructure
  std::vector<PolyBezierGroup> ComposeMicrostructureAndDerivatives() const {
    // Loop over the externnal splines in the group
    const PolyBezierGroup deformation_function =
        deformation_function_generator.Create();

    // Write it out (if required)
    utils::Export::GuessByExtension(deformation_function,
                                    "deformation_function_generator.xml");

    // Initialize return value
    std::vector<PolyBezierGroup> return_value(
        ValueField::kNumberOfSuperParameters + 1u,
        PolyBezierGroup(Microtile::kNumberOfSplines *
                        deformation_function.size()));

    // Every spline in the mirostructure defines a specific range of the
    // total spline group. We consider an even splitting, meaning that that
    // the parametric sections have the same size. As a first step, the
    // current field needs to be identified.
    const double dx =
        1. / static_cast<double>(
                 deformation_function_generator.GetNumberOfSegments()[0]);
    const double dy =
        1. / static_cast<double>(
                 deformation_function_generator.GetNumberOfSegments()[1]);
    const double dz =
        1. / static_cast<double>(
                 deformation_function_generator.GetNumberOfSegments()[2]);

#ifdef ENABLE_OPEN_MP_PARALLEL
#pragma omp parallel for
#endif
    for (std::size_t i_def_x = 0;
         i_def_x < deformation_function_generator.GetNumberOfSegments()[0];
         i_def_x++) {
      for (std::size_t i_def_y{};
           i_def_y < deformation_function_generator.GetNumberOfSegments()[1];
           i_def_y++) {
        for (std::size_t i_def_z{};
             i_def_z < deformation_function_generator.GetNumberOfSegments()[2];
             i_def_z++) {
          std::vector<PolyBezierGroup> microtile_vector;
          if constexpr (Microtile::HAS_CLOSING_TILE_DEFINITION) {
            // If a closing tile is defined we can create watertight structures
            if (i_def_y == 0) {
              const auto transformed_points = TransformPoints_(
                  Point3D{i_def_x * dx, i_def_y * dy, i_def_z * dz},
                  Point3D{(i_def_x + 1) * dx, (i_def_y + 1) * dy,
                          (i_def_z + 1) * dz},
                  Microtile::kClosingTileEvaluationPoints);
              // Evaluate microtile based on the generator provided by
              microtile_vector = Microtile::GenerateMicrostructureClosingTile(
                  Microtile::ClosingFace::X_MAX,
                  value_field.Evaluate(transformed_points));
            } else if (i_def_x ==
                       deformation_function_generator.GetNumberOfSegments()[1] -
                           1) {
              const auto transformed_points = TransformPoints_(
                  Point3D{i_def_x * dx, i_def_y * dy, i_def_z * dz},
                  Point3D{(i_def_x + 1) * dx, (i_def_y + 1) * dy,
                          (i_def_z + 1) * dz},
                  Microtile::kClosingTileEvaluationPoints);
              // Evaluate microtile based on the generator provided by
              microtile_vector = Microtile::GenerateMicrostructureClosingTile(
                  Microtile::ClosingFace::X_MAX,
                  value_field.Evaluate(transformed_points));
            } else {
              // If no closing tile is defined, then we set normal tiles in the
              // boundary elements
              const auto transformed_points = TransformPoints_(
                  Point3D{i_def_x * dx, i_def_y * dy, i_def_z * dz},
                  Point3D{(i_def_x + 1) * dx, (i_def_y + 1) * dy,
                          (i_def_z + 1) * dz},
                  Microtile::kEvaluationPoints);
              // Evaluate microtile based on the generator provided by
              microtile_vector = Microtile::GenerateMicrostructureDerivatives(
                  value_field.Evaluate(transformed_points));
            }
          } else {
            // If no closing tile is defined, then we set normal tiles in the
            // boundary elements
            const auto transformed_points = TransformPoints_(
                Point3D{i_def_x * dx, i_def_y * dy, i_def_z * dz},
                Point3D{(i_def_x + 1) * dx, (i_def_y + 1) * dy,
                        (i_def_z + 1) * dz},
                Microtile::kEvaluationPoints);
            // Evaluate microtile based on the generator provided by Microtile
            microtile_vector = Microtile::GenerateMicrostructureDerivatives(
                value_field.Evaluate(transformed_points));
          }

          // Compose the microstructure in first position of vector
          const std::size_t combined_index =
              i_def_x +
              deformation_function_generator.GetNumberOfSegments()[0] *
                  i_def_y +
              deformation_function_generator.GetNumberOfSegments()[0] *
                  deformation_function_generator.GetNumberOfSegments()[1] *
                  i_def_z;

          // Construct the actual Composition (microstructure)
          const auto composed_microstructure_tile =
              deformation_function[combined_index].Compose(microtile_vector[0]);
          for (unsigned i_spline_in_tile{};
               i_spline_in_tile < Microtile::kNumberOfSplines;
               i_spline_in_tile++) {
            return_value[0][combined_index * Microtile::kNumberOfSplines +
                            i_spline_in_tile] =
                composed_microstructure_tile[i_spline_in_tile];
          }

          /// Construct the corresponding derivatives
          // Precompute the derivatives of the deformation function as they are
          // required multiple times (for every internal control points of the
          // value field funciton)
          // Here it is hard-coded for the two dimensional case
          auto microstructure_derivative_first_par_dim =
              deformation_function[combined_index]
                  .DerivativeWRTParametricDimension(0)
                  .Compose(microtile_vector[0]);
          auto microstructure_derivative_second_par_dim =
              deformation_function[combined_index]
                  .DerivativeWRTParametricDimension(1)
                  .Compose(microtile_vector[0]);
          auto microstructure_derivative_third_par_dim =
              deformation_function[combined_index]
                  .DerivativeWRTParametricDimension(2)
                  .Compose(microtile_vector[0]);
          for (std::size_t i_derivative{};
               i_derivative < value_field.kNumberOfSuperParameters;
               i_derivative++) {
            const auto microstructure_derivative =
                (microstructure_derivative_first_par_dim.MultiplyComponentwise(
                     microtile_vector[i_derivative + 1].ExtractDimension(0)))
                    .AddComponentwise(microstructure_derivative_second_par_dim
                                          .MultiplyComponentwise(
                                              microtile_vector[i_derivative + 1]
                                                  .ExtractDimension(1)))
                    .AddComponentwise(microstructure_derivative_third_par_dim
                                          .MultiplyComponentwise(
                                              microtile_vector[i_derivative + 1]
                                                  .ExtractDimension(2)));
            // Assign to microstructure
            for (unsigned i_spline_in_tile{};
                 i_spline_in_tile < Microtile::kNumberOfSplines;
                 i_spline_in_tile++) {
              return_value[i_derivative + 1]
                          [combined_index * Microtile::kNumberOfSplines +
                           i_spline_in_tile] =
                              microstructure_derivative[i_spline_in_tile];
            }
          }
        }
      }
    }
    return return_value;
  }

  // Compose Derivatives
};

#endif  // EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_MICROSTRUCTURE_GENERATOR_HPP
