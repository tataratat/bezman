#ifndef SRC_UTILS_BASE64_HPP
#define SRC_UTILS_BASE64_HPP

#include <array>
#include <vector>

#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/type_traits/is_point.hpp"

namespace beziermanipulation::utils {

/**
 * @brief Encode for base64 export
 *
 * Python readable base64 export for double values. This allows to store doubles
 * exactly in textformat
 *
 */
class Base64 {
 private:
  /// Alias for one byte type
  using ByteRepresentation = unsigned char;

  /// Look up table
  static constexpr std::array<char, 64> char_encode_table{
      'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
      'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
      'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
      'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
      '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'};

 public:
  /**
   * @brief Actual encoding routine
   *
   * @tparam BaseType type of individual data entries
   * @param data_vector data to be encoded
   * @return std::string encoded data
   */
  template <typename BaseType, std::enable_if_t<!type_traits::isPoint_v<BaseType>, int> = 0>
  static std::string Encode(const std::vector<BaseType>& data_vector) {
    const ByteRepresentation* vector_as_bytes =
        reinterpret_cast<const ByteRepresentation*>(&data_vector[0]);

    // Number of bytes for an entry
    const std::size_t length_of_entry{sizeof(BaseType{})};
    // Minimum number of bytes required
    const std::size_t minimum_n_bytes_required =
        length_of_entry * data_vector.size();
    // Number of bytes must be divisible by three
    const std::size_t additional_padding_bytes =
        (3 - minimum_n_bytes_required % 3) % 3;
    // Required groups of three
    const std::size_t number_of_groups =
        (minimum_n_bytes_required + additional_padding_bytes) / 3;

    // Initialize return value
    std::string encoded_string;
    encoded_string.resize(number_of_groups * 4);

    // Loop over bytes and decode them
    for (std::size_t i_group{}; i_group < number_of_groups; i_group++) {
      const std::size_t buffer_index = i_group * 3;
      std::array<ByteRepresentation, 3> buffer{};
      buffer[0] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 0]
                      : 0;
      buffer[1] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 1]
                      : 0;
      buffer[2] = buffer_index < minimum_n_bytes_required
                      ? vector_as_bytes[buffer_index + 2]
                      : 0;

      encoded_string[i_group * 4 + 0] =
          char_encode_table[((buffer[0] & 0xfc) >> 2)];
      encoded_string[i_group * 4 + 1] =
          char_encode_table[((buffer[0] & 0x03) << 4) +
                            ((buffer[1] & 0xf0) >> 4)];
      encoded_string[i_group * 4 + 2] =
          char_encode_table[((buffer[1] & 0x0f) << 2) +
                            ((buffer[2] & 0xc0) >> 6)];
      encoded_string[i_group * 4 + 3] =
          char_encode_table[((buffer[2] & 0x3f) << 0)];
    }

    // Replace trailing invalid data with =
    for (size_t i = 0; i < additional_padding_bytes; ++i) {
      encoded_string[number_of_groups * 4 - i - 1] = '=';
    }

    return encoded_string;
  }

  /**
   * @brief Overload for Point type (handy with splines)
   *
   */
  template <std::size_t dimension, typename Scalar>
  static std::string Encode(
      const std::vector<beziermanipulation::Point<dimension, Scalar>>&
          data_vector) {
    std::vector<Scalar> ctps_converted(dimension * data_vector.size());
    for (std::size_t i_point{}; i_point < data_vector.size(); i_point++) {
      for (std::size_t i_dim{}; i_dim < dimension; i_dim++) {
        ctps_converted[i_point * dimension + i_dim] =
            data_vector[i_point][i_dim];
      }
    }
    return Encode(ctps_converted);
  }
};
}  // namespace beziermanipulation::utils

#endif  // SRC_UTILS_BASE64_HPP
