// #include <random>

// Load tensorSuite
#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

int main() {
  // Define aliases
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

  // Define splines for composition
  std::vector<Point3D> surface_ctps{Point3D{0., 0., 0.}, Point3D{2., 1., 0.},
                                    Point3D{0., 2., 0.}, Point3D{2., 3., 0.}};
  std::vector<Point2D> line_ctps{Point2D{0., 0.}, Point2D{1., 0},
                                 Point2D{1., 1.}};
  std::array<std::size_t, 2> surface_degrees{1, 1};
  std::array<std::size_t, 1> line_degrees{2};
  BezierSpline<2, Point3D, double> surface =
      BezierSpline<2, Point3D, double>(surface_degrees, surface_ctps);
  BezierSpline<1, Point2D, double> line =
      BezierSpline<1, Point2D, double>(line_degrees, line_ctps);

  // Perform Composition
  const auto composed_spline = surface.compose(line);

  // Print Composition
  for (unsigned int i{}; i < composed_spline.NumberOfControlPoints; i++){
    std::cout << composed_spline.control_points[i] << std::endl;
  }

  composed_spline.evaluate(0.5);

  return 0;
}
