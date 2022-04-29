#ifndef UTILS_ALGORITHMS_SORT_HPP
#define UTILS_ALGORITHMS_SORT_HPP

#include <algorithm>
#include <numeric>
#include <vector>

namespace beziermanipulation::utils::algorithms {

/*
 * Sort Vector using lambda expressions
 * https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
 */
template <typename T>
std::vector<std::size_t> IndexListSort(const std::vector<T>& v) {
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&v](size_t i1, size_t i2) { return v[i1] < v[i2]; });

  return idx;
}

}  // namespace beziermanipulation::utils::algorithms
#endif  // UTILS_ALGORITHMS_SORT_HPP
