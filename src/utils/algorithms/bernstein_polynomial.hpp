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

#ifndef UTILS_ALGORITHMS_BERNSTEIN_POLYNOMIAL_HPP
#define UTILS_ALGORITHMS_BERNSTEIN_POLYNOMIAL_HPP

#include <cassert>
#include <vector>

#include "bezman/src/utils/fastbinomialcoefficient.hpp"

namespace bezman::utils::algorithms {
/**
 * @brief
 *
 * @tparam ScalarType
 */
template <typename ScalarType = double>
class BernsteinPolynomial {
 public:
  /**
   * @brief Evaluate Bernstein polynomials
   *
   * Efficient calculations of Bernstein polynomials to avoid power expressions
   *
   * @param value   position to be evaluated
   * @param degree  polynomial degrees
   * @return constexpr std::vector<ScalarType> vector of evaluations of the
   * different basis functions
   */
  static inline constexpr std::vector<ScalarType> Evaluate(
      const std::size_t& degree, const ScalarType& value) {
    // Check if degree is available
    assert(degree < MAX_BINOMIAL_DEGREE);

    // Init Return Type
    std::vector<ScalarType> evaluations;

    // Preallocate space to vectors
    evaluations.reserve(degree + 1);

    // Loop over parametric dimensions to precalculate basisfunctions
    ScalarType cfactor = static_cast<ScalarType>(1.);
    const ScalarType inv_value = cfactor - value;

    // Loop over dimensionwise ctps
    for (std::size_t i_basis{0}; i_basis < degree + 1; i_basis++) {
      evaluations.push_back(
          static_cast<ScalarType>(
              utils::FastBinomialCoefficient::choose(degree, i_basis)) *
          cfactor);
      cfactor *= value;
    }
    // Backwards assignment
    cfactor = static_cast<ScalarType>(1.);
    for (std::size_t i_basis{0}; i_basis < degree + 1; i_basis++) {
      evaluations[degree - i_basis] *= cfactor;
      cfactor *= inv_value;
    }
    return evaluations;
  }

  /**
   * @brief
   *
   * @param value
   * @param degree
   * @param deriv
   * @return constexpr std::vector<ScalarType>
   */
  static inline constexpr std::vector<ScalarType> EvaluateDerivative(
      const std::size_t& degree, const ScalarType& value,
      const std::size_t& deriv) {
    // Check if degree is available
    assert(degree < MAX_BINOMIAL_DEGREE);

    // Init Return Type
    std::vector<ScalarType> evaluations;

    // Determine derivatives
    if (degree < deriv) {
      evaluations.resize(degree + 1);
      return evaluations;
    } else if (deriv == 0) {
      return Evaluate(degree, value);
    } else {
      // based on doi:10.1155/2011/829543, section 3
      // Retrive low order basis functions
      const std::vector<ScalarType> basis = Evaluate(degree - deriv, value);
      // Auxiliary value
      const ScalarType aux_factor = [&degree, &deriv]() {
        ScalarType value{1};
        for (std::size_t i_deg{degree}; i_deg > (degree - deriv); i_deg--) {
          value *= i_deg;
        }
        return value;
      }();
      // Start actual calculation using the papers notation
      evaluations.resize(degree + 1);
      for (std::size_t i_basis{0}; i_basis < degree + 1; i_basis++) {
        const std::size_t k_min{static_cast<std::size_t>(
            std::max(0, static_cast<int>(i_basis + deriv - degree)))};
        for (std::size_t k{k_min}; k <= std::min(i_basis, deriv); k++) {
          // Check if is even using binary operation
          const ScalarType eval_v =
              utils::FastBinomialCoefficient::choose(deriv, k) *
              basis[i_basis - k];
          evaluations[i_basis] += ((k + deriv) & 1) ? -eval_v : eval_v;
        }
        evaluations[i_basis] *= aux_factor;
      }
      return evaluations;
    }
  }
};
}  // namespace bezman::utils::algorithms

#endif  // UTILS_ALGORITHMS_BERNSTEIN_POLYNOMIAL_HPP
