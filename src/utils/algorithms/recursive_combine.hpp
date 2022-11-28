#ifndef UTILS_ALGORITHMS_RECURSIVE_COMBINE_HPP
#define UTILS_ALGORITHMS_RECURSIVE_COMBINE_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <vector>

namespace bezman::utils::algorithms {

/**
 * @brief Actual Implementation of Recursive Combine (should not be used)
 */
template <std::size_t kDepth, std::size_t parametric_dimension,
          typename ValueType>
constexpr inline void RecursiveCombine_(
    const std::array<std::vector<ValueType>, parametric_dimension>& factors,
    std::vector<ValueType>& result, const ValueType& c_value) {
  static_assert(kDepth < parametric_dimension,
                "Implementation error, recursion loop to deep!");
  for (std::size_t i_entry{}; i_entry < factors[kDepth].size(); i_entry++) {
    if constexpr (kDepth == 0) {
      result.push_back(c_value * factors[kDepth][i_entry]);
    } else {
      RecursiveCombine_<static_cast<std::size_t>(kDepth - 1)>(
          factors, result, c_value * factors[kDepth][i_entry]);
    }
  }
}

/**
 * @brief Using efficient recursion at compile time we can avoid multple
 * calculations of the same product
 *
 * This is done in the following way:
 * For a three dimensional p
 */
template <std::size_t parametric_dimension, typename ValueType>
constexpr std::vector<ValueType> RecursiveCombine(
    const std::array<std::vector<ValueType>, parametric_dimension>& factors) {
  // Precalculate required entries
  std::size_t n_entries{1};
  std::for_each(
      factors.begin(), factors.end(),
      [&](const std::vector<ValueType>& ii) { n_entries *= ii.size(); });

  // Init return type and reserve memory
  std::vector<ValueType> result{};
  result.reserve(n_entries);

  // Start computation
  RecursiveCombine_<parametric_dimension - 1uL>(factors, result,
                                                static_cast<ValueType>(1));

  return result;
}

}  // namespace bezman::utils::algorithms

#endif  // UTILS_ALGORITHMS_RECURSIVE_COMBINE_HPP
