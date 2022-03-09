// #include <random>

// Load tensorSuite
#include "bezierManipulation/src/bezier_spline.hpp"

using namespace beziermanipulation;

int main() {
  // Define aliases
  using Point3D = Point<3, double>;
  using Point2D = Point<2, double>;

  // Define splines for composition
  std::vector<Point3D> surface_ctps{Point3D{1., 1., 0.}, Point3D{3., 2., 0.},
                                    Point3D{1., 3., 0.}, Point3D{3., 4., 0.}};
  std::vector<Point2D> line_ctps{Point2D{0., 0.}, Point2D{1., 0},
                                 Point2D{1., 1.}};
  std::vector<Point2D> inner_surface_ctps{Point2D{0.2, 0.}, Point2D{0.8, 0},
                                          Point2D{0.2, 0.8}, Point2D{0.2, 0.8}};
  std::array<std::size_t, 2> surface_degrees{1, 1};
  std::array<std::size_t, 1> line_degrees{2};
  BezierSpline<2, Point3D, double> surface =
      BezierSpline<2, Point3D, double>(surface_degrees, surface_ctps);
  BezierSpline<1, Point2D, double> line =
      BezierSpline<1, Point2D, double>(line_degrees, line_ctps);
  BezierSpline<2, Point2D, double> inner_surface =
      BezierSpline<2, Point2D, double>(surface_degrees, inner_surface_ctps);

  // Perform Composition
  const auto composed_line_spline = surface.compose(line);
  const auto composed_surface_spline = surface.compose(inner_surface);

  // Print Composition
  std::cout << "Surface-Line-Composition\n";
  for (unsigned int i{}; i < composed_line_spline.NumberOfControlPoints; i++) {
    std::cout << composed_line_spline.control_points[i] << std::endl;
  }

  std::cout << "Surface-Surface-Composition\n";
  for (unsigned int i{}; i < composed_surface_spline.NumberOfControlPoints;
       i++) {
    std::cout << composed_surface_spline.control_points[i] << std::endl;
  }

  // Evaluate some points for testing
  for (int i{}; i < 10; i++) {
    const double x{static_cast<double>(rand()) / static_cast<double>(RAND_MAX)};
    const auto composed_spline_position = composed_line_spline.evaluate(x);
    const auto composed_spline_expectation = surface.evaluate(line.evaluate(x));
    std::cout << "Composed spline evaluation at : " << x << "\t yields "
              << composed_spline_position ;
    std::cout << "\tExpected position is : " << composed_spline_expectation << "\n";
  }
    composed_line_spline.evaluate(0.5);
    return 0;
  }
