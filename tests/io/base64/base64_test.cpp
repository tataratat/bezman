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

#include "bezman/src/utils/base64.hpp"

#include <gtest/gtest.h>

#include <array>
#include <iostream>
#include <vector>

#include "bezman/src/point.hpp"

using namespace bezman::utils;

namespace bezman::tests::io::base64 {

class Base64ExportSuite : public ::testing::Test {
 public:
  std::vector<int> int_vect_one{1, 4, -3, 5, 10, 2};
  std::vector<int> int_vect_two{1, 4, -3, 5};
  std::vector<float> float_vect_one{3.0, 4., 0.2, -.4};
  std::vector<double> double_vect_one{3.0, 4., 0.2, -.4};

  std::vector<bezman::Point<3, double>> point_vector{
      bezman::Point<3, double>{1., 2., 1.},
      bezman::Point<3, double>{-0.65, 0., 9.}};
};

/*
 * Encode and Decode some vectors
 */
TEST_F(Base64ExportSuite, ExportImportTest) {
  // Check first against python export in b64
  EXPECT_EQ(std::string("AQAAAAQAAAD9////BQAAAAoAAAACAAAA"),
            Base64::Encode(int_vect_one));
  // Check if encoding and decoding results in the same thing
  EXPECT_EQ(int_vect_one, Base64::Decode<int>(Base64::Encode(int_vect_one)));
  EXPECT_EQ(int_vect_two, Base64::Decode<int>(Base64::Encode(int_vect_two)));
  EXPECT_EQ(float_vect_one,
            Base64::Decode<float>(Base64::Encode(float_vect_one)));
  EXPECT_EQ(double_vect_one,
            Base64::Decode<double>(Base64::Encode(double_vect_one)));
}

TEST_F(Base64ExportSuite, ExportImportTestPoints) {
  // Points
  const auto test_point_ctps =
      Base64::Decode<bezman::Point<3ul, double>>(Base64::Encode(point_vector));
  EXPECT_EQ(point_vector, test_point_ctps);
}

}  // namespace bezman::tests::io::base64