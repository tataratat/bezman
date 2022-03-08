#ifndef SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
#define SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP

#include <array>

#include "bezierManipulation/src/utils/binomialcoefficientlookuptablecreator.hpp"

#ifndef MAX_BINOMIAL_DEGREE
#define MAX_BINOMIAL_DEGREE 30u
#endif


namespace beziermanipulation::utils{

// Global Lookup
class FastBinomialCoefficient {
 private:
  using ScalarType = std::size_t;
  constexpr static auto n_values =
      BinomialCoefficientLookupTableCreator<ScalarType,
                                            MAX_BINOMIAL_DEGREE>::n_values;
  constexpr static std::array<ScalarType, n_values> look_up =
      BinomialCoefficientLookupTableCreator<
          ScalarType, MAX_BINOMIAL_DEGREE>::calculate_table();

 public:
  constexpr static ScalarType choose(const unsigned int n,
                                     const unsigned int i) {
    assert(n < MAX_BINOMIAL_DEGREE);
    return look_up[n * (n + 1) / 2 + i];
  }
};
}  // namespace beziermanipulation::utils

#endif  // SRC_UTILS_FASTBINOMIALCOEFFICIENT_HPP
