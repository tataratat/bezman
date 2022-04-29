#ifndef UTILS_ALGORITHMS_INT_POWER_HPP
#define UTILS_ALGORITHMS_INT_POWER_HPP

namespace beziermanipulation::utils::algorithms {

/**
 * @brief Recursive function of integer powers
 */
template <typename IntegralType,
          typename E = std::enable_if_t<std::is_integral_v<IntegralType>>>
constexpr IntegralType IntPower(const IntegralType& base,
                                const IntegralType& power) {
  if (power == 0) return 1;
  if (power == 1) return base;

  int tmp = IntPower(base, power / 2);
  if ((power % 2) == 0) {
    return tmp * tmp;
  } else {
    return base * tmp * tmp;
  }
}
}  // namespace beziermanipulation::utils::sort

#endif  // UTILS_ALGORITHMS_INT_POWER_HPP
