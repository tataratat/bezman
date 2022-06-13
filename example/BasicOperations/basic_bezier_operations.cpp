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

// #include <random>

// Load tensorSuite
#include "bezman/src/bezier_spline.hpp"

using namespace bezman;

int main() {
  std::cout << "Working on the Bezier Splines\n";
  const unsigned int n_x{1}, n_y{1}, n_z{1};

  auto spline_test = BezierSpline<3u, Point<3, double>, double>(
      std::array<std::size_t, 3>{n_x, n_y, n_z});

  for (unsigned int i{0}; i <= n_x; i++) {
    for (unsigned int j{0}; j <= n_y; j++) {
      for (unsigned int k{0}; k <= n_z; k++) {
        spline_test.ControlPoint(i, j, k) =
            Point<3, double>((double)i, (double)j, (double)k);
      }
    }
  }

  BezierSpline<3u, Point<3, double>, double> spline_test_copy{spline_test};

  const auto sum_spline = spline_test * spline_test_copy;

  for (unsigned int i{0}; i < sum_spline.NumberOfControlPoints; i++) {
    std::cout << sum_spline.control_points[i];
    std::cout << " : " << i << std::endl;
  }

  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        y{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)},
        z{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto eval_point = sum_spline.Evaluate(Point<3, double>{x, y, z});
    const auto eval_point1 = spline_test.Evaluate(Point<3, double>{x, y, z});
    const auto eval_point2 =
        spline_test_copy.Evaluate(Point<3, double>{x, y, z});

    std::cout << std::setw(5) << std::setprecision(3) << "Result :\t"
              << eval_point << "\tVec1 :\t" << eval_point1 << "\tVec2 :\t"
              << eval_point2 << "\tExpected Result :\t"
              << eval_point1 * eval_point2 << std::endl;
  }

  // // Second Test
  using Point3D = Point<3, double>;

  // Provide some data to form splines.
  std::vector<Point3D> line1ctps{Point3D{0., 0., 0.}, Point3D{0., 0.5, 0.},
                                 Point3D{0., 1., 0.}};
  std::vector<Point3D> line2ctps{Point3D{0., 0., 0.}, Point3D{0.5, 0., 0.},
                                 Point3D{1., 0., 0.}};
  std::vector<Point3D> line3ctps{Point3D{0., 0., 0.}, Point3D{0.5, 0.5, 0.},
                                 Point3D{1., 1., 0.}};
  std::array<std::size_t, 1> degrees{2};

  // Some Lines
  BezierSpline<1, Point3D, double> line1 =
      BezierSpline<1, Point3D, double>(degrees, line1ctps);
  BezierSpline<1, Point3D, double> line2 =
      BezierSpline<1, Point3D, double>(degrees, line2ctps);
  BezierSpline<1, Point3D, double> line3 =
      BezierSpline<1, Point3D, double>(degrees, line3ctps);
  const auto test0 =
      (line1.OrderElevateAlongParametricDimension(0) + line2);
  const auto test1 =
      line3.OrderElevateAlongParametricDimension(0);
  std::cout
      << ((line1.OrderElevateAlongParametricDimension(0) + line2) ==
          line3.OrderElevateAlongParametricDimension(0));
  return 0;
}
