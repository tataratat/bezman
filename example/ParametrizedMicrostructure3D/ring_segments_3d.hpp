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

#ifndef EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_RING_SEGMENTS_3D_HPP
#define EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_RING_SEGMENTS_3D_HPP

/**
 * @brief Deformation function defining the outer contour of a spline
 *
 * Approximates a circle with a given number of segments in each parametric
 * dimension. Circles are approximated with second order B-Spline segments.
 *
 */
class RingSegments3D {
 private:
  // Aliases
  using ADT = utils::computational_differentiation::AlgoDiffType<double>;
  using PointADT3D = Point<3, ADT>;
  using Point3D = Point<3, double>;
  using Bezier = BezierSpline<3, Point3D, double>;
  using PolyBezierGroup = BezierGroup<Bezier>;

  /// Precalculated values
  const double PI = std::acos(-1.);

  // Default values are set here
  double innerR{1.}, outerR{2.}, arc_degrees{PI * 0.25}, depth{2.};
  std::array<std::size_t, 3> numberOfSegments{1, 1, 1};

 public:
  // All Setter Methods return a reference to the object so that multiple
  // values can be changed in a single line for convenience Quarter Circle
  // Dimensions

  /// Set inner radius
  constexpr RingSegments3D& SetInnerRadius(const double r) {
    innerR = r;
    return (*this);
  }

  /// Set outer radius
  constexpr RingSegments3D& SetOuterRadius(const double r) {
    outerR = r;
    return (*this);
  }

  /// Set Radius
  constexpr RingSegments3D& SetArcLength(const double r) {
    arc_degrees = r;
    return (*this);
  }

  /// Set Radius
  constexpr RingSegments3D& SetDepth(const double d) {
    depth = d;
    return (*this);
  }

  /// Set number of segments
  constexpr RingSegments3D& SetNumberOfSegments(
      const std::array<std::size_t, 3>& number_of_segments) {
    numberOfSegments = number_of_segments;
    return (*this);
  }

  /// Get number of segments
  constexpr const std::array<std::size_t, 3>& GetNumberOfSegments() const {
    return numberOfSegments;
  }

  /// Default Constructor
  RingSegments3D() = default;

  /**
   * @brief Create Circle Deformation function with given number of segments
   * per parametric dimension
   */
  PolyBezierGroup Create() const {
    // split planes to segment circle in radial direction (even samples)
    std::vector<double> y_knot_lines(numberOfSegments[1] - 1);
    const double y_seg_length = 1. / static_cast<double>(numberOfSegments[1]);
    for (unsigned int i_y_segment{1}; i_y_segment < numberOfSegments[1];
         i_y_segment++) {
      y_knot_lines[i_y_segment - 1] = y_seg_length * i_y_segment;
    }

    // Initialize return value
    PolyBezierGroup ringsegments{numberOfSegments[0] * numberOfSegments[1] *
                                 numberOfSegments[2]};

    // Precompute values that are required multiple times
    const std::array<std::size_t, 3> degrees{2, 1, 1};
    const double degrees_per_segment = arc_degrees / numberOfSegments[0];
    const double depth_per_segment = depth / numberOfSegments[2];

    const double excentricity_of_middle_points =
        1. / std::sin(PI / 2. - degrees_per_segment / 2.);
    for (std::size_t i_z_segment{}; i_z_segment < numberOfSegments[2];
         i_z_segment++) {
      // Precompute depths
      const double z_start = depth_per_segment * i_z_segment;
      const double z_end = depth_per_segment * (i_z_segment + 1);
      for (std::size_t i_x_segment{}; i_x_segment < numberOfSegments[0];
           i_x_segment++) {
        const double startsin = std::sin(degrees_per_segment * i_x_segment);
        const double startcos = std::cos(degrees_per_segment * i_x_segment);
        const double middlesin =
            std::sin(degrees_per_segment * (i_x_segment + .5));
        const double middlecos =
            std::cos(degrees_per_segment * (i_x_segment + .5));
        const double endsin =
            std::sin(degrees_per_segment * (i_x_segment + 1.));
        const double endcos =
            std::cos(degrees_per_segment * (i_x_segment + 1.));

        // Define the control points
        std::vector<Point3D> ctps{
            Point3D{startcos * outerR, startsin * outerR, z_start},
            Point3D{middlecos * excentricity_of_middle_points * outerR,
                    middlesin * excentricity_of_middle_points * outerR,
                    z_start},
            Point3D{endcos * outerR, endsin * outerR, z_start},
            Point3D{startcos * innerR, startsin * innerR, z_start},
            Point3D{middlecos * excentricity_of_middle_points * innerR,
                    middlesin * excentricity_of_middle_points * innerR,
                    z_start},
            Point3D{endcos * innerR, endsin * innerR, z_start},
            Point3D{startcos * outerR, startsin * outerR, z_end},
            Point3D{middlecos * excentricity_of_middle_points * outerR,
                    middlesin * excentricity_of_middle_points * outerR, z_end},
            Point3D{endcos * outerR, endsin * outerR, z_end},
            Point3D{startcos * innerR, startsin * innerR, z_end},
            Point3D{middlecos * excentricity_of_middle_points * innerR,
                    middlesin * excentricity_of_middle_points * innerR, z_end},
            Point3D{endcos * innerR, endsin * innerR, z_end}};
        // Define and split up circle segment
        const auto CircleWedge = Bezier{degrees, ctps}
                                     //.OrderElevateAlongParametricDimension(1)
                                     //.OrderElevateAlongParametricDimension(2)
                                     .SplitAtPosition(y_knot_lines, 1);
        // Assign splines to Splinegroup
        const int index_offset = i_x_segment + numberOfSegments[0] *
                                                   numberOfSegments[1] *
                                                   i_z_segment;
        for (std::size_t i_y_segment{}; i_y_segment < numberOfSegments[1];
             i_y_segment++) {
          ringsegments[index_offset + numberOfSegments[0] * i_y_segment] =
              CircleWedge[i_y_segment];
        }
      }
    }
    return ringsegments;
  }
};

#endif  // EXAMPLE_PARAMETRIZED_MICROSTRUCTURE3D_RING_SEGMENTS_3D_HPP
