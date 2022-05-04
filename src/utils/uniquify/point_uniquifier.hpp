#ifndef UTILS_UNIQUIFY_POINT_UNIQUIFIER_HPP
#define UTILS_UNIQUIFY_POINT_UNIQUIFIER_HPP

#include <array>
#include <cassert>
#include <vector>

#include "bezierManipulation/src/bezier_spline_group.hpp"
#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/algorithms/hypercube.hpp"
#include "bezierManipulation/src/utils/algorithms/sort.hpp"

namespace beziermanipulation::utils::uniquify {

/**
 * @brief Takes a set of points and determines point-connectivity using a metric
 *
 * Duplicate Points are not eliminated, assuming that aa maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam physical_dimension   Entries in Point
 * @tparam ScalarType           Type of Entry in Point
 * @param face_center_points    List of intersection points (Mid-Point of
 *                              face-vertices)
 * @param orientation_metric    Vector along which the points are ordered
 *                              prior to their uniquification
 * @param tolerance             Tolerance for vector contraction
 * @return Connectivity         Element Connectivity
 */
template <std::size_t physical_dimension, typename ScalarType,
          std::size_t number_of_element_faces>
auto FindConnectivity(
    const std::vector<beziermanipulation::Point<
        physical_dimension, ScalarType>>& face_center_points,
    const beziermanipulation::Point<physical_dimension, ScalarType>
        orientation_metric,
    const std::array<std::size_t, number_of_element_faces>& opposite_face_list,
    const ScalarType tolerance = 1e-5) {
  // Check if number of faces is a divisor of the point list length
  assert(face_center_points.size() % number_of_element_faces == 0);

  // Assure Metric is normed and non-zero
  assert(orientation_metric.EuclidianNorm() > 0);
  const beziermanipulation::Point<physical_dimension, ScalarType>
      normed_orientation_metric =
          orientation_metric *
          (static_cast<ScalarType>(1.) / orientation_metric.EuclidianNorm());

  // Store information in Auxiliary Values
  const std::size_t n_total_points{face_center_points.size()};
  const ScalarType tolerance_squared{tolerance * tolerance};

  // Init connectivity and metric value
  // (-1 : boundary, -2 : untouched)
  // {in c++20 this expression could be constexpr}
  const std::size_t number_of_elements =
      face_center_points.size() / number_of_element_faces;
  std::vector<std::array<int, number_of_element_faces>> connectivity(
      // Size of vector
      number_of_elements,
      // Lambda function to initialize an array with constant value in size of
      // face-number
      []() {
        std::array<int, number_of_element_faces> a{};
        a.fill(-2);
        return a;
      }());
  std::vector<ScalarType> scalar_metric(n_total_points);

  // Check Metric Dimension and Vector Size
  for (unsigned int i{}; i < n_total_points; i++) {
    scalar_metric[i] = normed_orientation_metric * face_center_points[i];
  }

  // Sort Metric Vector
  const auto metric_order_indices = algorithms::IndexListSort(scalar_metric);

  // Loop over points
  for (unsigned int lower_limit = 0; lower_limit < n_total_points - 1;
       lower_limit++) {
    // Loop over all points regardless of whether they have been touched or not,
    // and then check the validity of the connection Point already processed
    bool found_duplicate = false;
    // Now check allowed range for duplicates
    unsigned int upper_limit = lower_limit + 1;
    while (scalar_metric[metric_order_indices[upper_limit]] -
                   scalar_metric[metric_order_indices[lower_limit]] <
               tolerance &&
           upper_limit < n_total_points) {
      // Check if the two points are duplicates
      found_duplicate = (face_center_points[metric_order_indices[lower_limit]] -
                         face_center_points[metric_order_indices[upper_limit]])
                            .SquaredEuclidianNorm() < tolerance_squared;
      if (found_duplicate) {
        break;
      } else {
        upper_limit++;
      }
    }

    // Now we have to check if the connection is valid
    // 1. If another connection is found, that means, that the point it connects
    //    to has a higher index in the metric tensor. If the current point does
    //    already have a neighbor, that means that more than one point connect
    //    -> Error
    // 2. The point it connects to must be on opposite sides on the neighboring
    //    element, else there is an orientation problem in the mesh and the mesh
    //    can not be exported in mfem format, this check needs to be disabled
    //    for general connectivity
    const std::size_t id_start_point = metric_order_indices[lower_limit];
    const std::size_t element_id_start =
        id_start_point / number_of_element_faces;
    const std::size_t element_face_id_start =
        id_start_point - element_id_start * number_of_element_faces;
    // is that id_start_point % number_of_element_faces?

    if (found_duplicate) {
      // Calculate indices
      const std::size_t id_end_point = metric_order_indices[upper_limit];
      const std::size_t element_id_end = id_end_point / number_of_element_faces;
      const std::size_t element_face_id_end =
          id_end_point - element_id_end * number_of_element_faces;
      // Check 1. (@todo EXCEPTION)
      assert(connectivity[element_id_start][element_face_id_start] == -2);
      // Check 2. (@todo EXCEPTION)
      assert(opposite_face_list[element_face_id_start] == element_face_id_end);
      assert(opposite_face_list[element_face_id_end] == element_face_id_start);
      // If both tests passed, update connectivity
      connectivity[element_id_start][element_face_id_start] = element_id_end;
      connectivity[element_id_end][element_face_id_end] = element_id_start;
    } else {
      // set Boundary-ID
      if (connectivity[element_id_start][element_face_id_start] == -2) {
        connectivity[element_id_start][element_face_id_start] = -1;
      }
    }
  }

  // Treat last remaining point in scalar metric vector
  const std::size_t last_id = metric_order_indices[n_total_points - 1];
  const std::size_t last_element = last_id / number_of_element_faces;
  const std::size_t last_face =
      last_id - last_element * number_of_element_faces;
  if (connectivity[last_element][last_face] == -2) {
    connectivity[last_element][last_face] = -1;
  }
  return connectivity;
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
auto GetConnectivityForSplineGroup(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group) {
  // Current implementation is only made for bi- and trivariates
  static_assert((
                    // parametric_dimension == 3 ||
                    parametric_dimension == 2),
                "High-Dimensional and Line Patches not supported");

  // Array that stores opposite faces
  constexpr auto opposite_faces =
      algorithms::HyperCube<parametric_dimension>::GetOppositeFaces();

  // Create Face-Center-Point Vector
  std::vector<PhysicalPointType> face_edges(spline_group.size() *
                                            opposite_faces.size());

  // Start Loop
  // (Instead of using the mean of the face vertices, using the sum)
  const std::size_t number_of_splines = spline_group.size();
  const std::size_t number_of_element_faces = opposite_faces.size();
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    const auto face_vertices =
        algorithms::HyperCube<parametric_dimension>::ControlPointIndicesOnFace(
            spline_group[i_spline].GetDegrees());
    for (std::size_t i_face{}; i_face < number_of_element_faces; i_face++) {
      face_edges[i_spline * number_of_element_faces + i_face] =
          spline_group[i_spline].control_points[face_vertices[i_face][0]];
      for (std::size_t i_point{1}; i_point < face_vertices[0].size();
           i_point++) {
        face_edges[i_spline * number_of_element_faces + i_face] +=
            spline_group[i_spline]
                .control_points[face_vertices[i_face][i_point]];
      }
    }
  }

  // Get Connectivity
  return FindConnectivity(
      face_edges,
      // Metric for internal ordering
      spline_group.MaximumCorner() - spline_group.MinimumCorner(),
      opposite_faces);
}

}  // namespace beziermanipulation::utils::uniquify
#endif  // UTILS_UNIQUIFY_POINT_UNIQUIFIER_HPP
