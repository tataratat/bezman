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

#ifndef UTILS_ALGORITHMS_POINT_UNIQUIFIER_HPP
#define UTILS_ALGORITHMS_POINT_UNIQUIFIER_HPP

#include <array>
#include <cassert>
#include <map>
#include <tuple>
#include <vector>

#include "bezman/src/bezier_group.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/utils/algorithms/hypercube.hpp"
#include "bezman/src/utils/algorithms/sort.hpp"

namespace bezman::utils::algorithms {

/**
 * @brief Determine the connectivity from center-vertices, assuming nothing of
 * the underlying grid
 *
 * Duplicate Points are not eliminated, assuming that a maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam PhysicalPointType Type of Point coordinates
 * @tparam ScalarType Type determining the precision
 * @tparam parametric_dimension dimension of the object (e.g. surface in 3D)
 * @tparam boolean check_orientation to check if neighboring elements match
 *                          structured grid
 * @param face_center_vertices vertices in the centers of spline-surfaces
 * @param metric used for preordering the vertices along a line
 * @param tolerance tolerance (distance between two vertices that are joined)

 * @return connectivity as a std::vector<std::array<...>>
 */
template <std::size_t parametric_dimension, bool check_orientation,
          typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::ScalarType_>
auto FindConnectivityFromCenters(
    const std::vector<PhysicalPointType>& face_center_vertices,
    const PhysicalPointType& metric, const ScalarType tolerance) {
  // -- Auxiliary data --
  constexpr std::size_t kParametricDimensions_ = parametric_dimension;
  constexpr std::size_t number_of_element_faces = parametric_dimension * 2;
  using ElementConnectivityInfoT =
      std::array<std::size_t, number_of_element_faces>;
  const std::size_t number_of_patches =
      face_center_vertices.size() / number_of_element_faces;
  const ScalarType tolerance_squared{tolerance * tolerance};

  // Consistency check
  if (!(face_center_vertices.size() % number_of_element_faces == 0)) {
    Logger::TerminatingError(
        "Number of corner vertices invalid. Must be a multiple of the number "
        "of vertices per patch");
  }

  // Check if number of faces is a divisor of the point list length
  Logger::Logging("Determining connectivity by analyzing face centers");
  // Check for Wrong number of faces and center points
  const std::size_t number_of_center_vertices{face_center_vertices.size()};
  assert(number_of_center_vertices % number_of_element_faces == 0);

  // Assure Metric is normed and non-zero
  const PhysicalPointType normed_metric = [](PhysicalPointType metric) {
    if (metric.EuclidianNorm() < 1e-20) {
      Logger::Warning(
          "Metric has no length. Chose non-zero "
          "metric for ordering points");
      Logger::Warning("Fall back to default metric, which is {1., 1., ...}");
      metric.fill(1.);
    }
    return metric * (static_cast<ScalarType>(1.) / metric.EuclidianNorm());
  }(metric);
  // Init connectivity and metric value
  // (-1 : boundary, -2 : untouched)
  std::vector<std::array<std::size_t, number_of_element_faces>> connectivity(
      // Size of vector
      number_of_patches,
      // Lambda function to initialize an array with constant value in size
      // of face-number
      [&]() {
        ElementConnectivityInfoT a{};
        a.fill(static_cast<std::size_t>(-2));
        return a;
      }());

  std::vector<ScalarType> scalar_metric{};
  scalar_metric.reserve(number_of_center_vertices);

  // Check Metric Dimension and Vector Size
  for (unsigned int i{}; i < number_of_center_vertices; i++) {
    scalar_metric.push_back(normed_metric * face_center_vertices[i]);
  }

  // Sort Metric Vector
  const auto metric_order_indices = algorithms::IndexListSort(scalar_metric);

  // Loop over points
  for (unsigned int lower_limit = 0;
       lower_limit < number_of_center_vertices - 1; lower_limit++) {
    // Loop over all points regardless of whether they have been touched or not,
    // and then check the validity of the connection Point already processed
    bool found_duplicate = false;
    // Now check allowed range for duplicates
    unsigned int upper_limit = lower_limit + 1;
    while (upper_limit < number_of_center_vertices &&
           scalar_metric[metric_order_indices[upper_limit]] -
                   scalar_metric[metric_order_indices[lower_limit]] <
               tolerance) {
      // Check if the two points are duplicates
      found_duplicate =
          (face_center_vertices[metric_order_indices[lower_limit]] -
           face_center_vertices[metric_order_indices[upper_limit]])
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
      if (connectivity[element_id_start][element_face_id_start] !=
          static_cast<std::size_t>(-2)) {
        Logger::TerminatingError(
            "Connectivity connection is invalid. "
            "Found conflicting interceptions");
      }

      // Check 2. (@todo EXCEPTION)
      // TODO check if mfem format is used for the output -> if not do not check
      if constexpr (check_orientation) {
        constexpr auto opposite_face_list =
            HyperCube<kParametricDimensions_>::GetOppositeFaces();
        if (opposite_face_list[element_face_id_start] != element_face_id_end) {
          Logger::TerminatingError("Orientation Problem for MFEM-mesh output.");
          // @todo In order to get the connectivity only, this check needs to be
          // performed, a boolean value needs to be switched, but the
          // connectivity is still returned
        }
#ifndef NDEBUG
        if (opposite_face_list[element_face_id_end] != element_face_id_start) {
          Logger::TerminatingError("Orientation Problem for MFEM-mesh output.");
        }
#endif
      }

      // If both tests passed, update connectivity
      connectivity[element_id_start][element_face_id_start] = element_id_end;
      connectivity[element_id_end][element_face_id_end] = element_id_start;
    } else {
      // set Boundary-ID
      if (connectivity[element_id_start][element_face_id_start] ==
          static_cast<std::size_t>(-2)) {
        connectivity[element_id_start][element_face_id_start] =
            static_cast<std::size_t>(-1);
      }
    }
  }

  // Treat last remaining point in scalar metric vector
  const std::size_t last_id =
      metric_order_indices[number_of_center_vertices - 1];
  const std::size_t last_element = last_id / number_of_element_faces;
  const std::size_t last_face =
      last_id - last_element * number_of_element_faces;
  if (connectivity[last_element][last_face] == static_cast<std::size_t>(-2)) {
    connectivity[last_element][last_face] = static_cast<std::size_t>(-1);
  }
  Logger::Logging("Found " + std::to_string(last_id) + " connections for " +
                  std::to_string(number_of_center_vertices) + " faces");
  return connectivity;
}

/**
 * @brief Determine the connectivity from corner-vertices, assuming nothing of
 * the underlying grid
 *
 * Duplicate Points are not eliminated, assuming that a maximum of two points
 * are equivalent. If this is not the case an exception is thrown. In theory
 * this has complexity O(nlogn) whereas a KDTree has complexity O(n (logn)^dim).
 *
 * @tparam PhysicalPointType Type of Point coordinates
 * @tparam ScalarType Type determining the precision
 * @tparam parametric_dimension dimension of the object (e.g. surface in 3D)
 * @param corner_vertices Corner Vertices that are extracted from the splines
 * @param metric used for preordering the vertices along a line
 * @param tolerance tolerance (distance between two vertices that are joined)
 * @param check_orientation boolean to check if neighboring elements match
 *                          structured grid
 * @return connectivity as a vector<array<...>>
 */
template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::ScalarType_>
auto FindConnectivityFromCorners(
    const std::vector<PhysicalPointType>& corner_vertices,
    const PhysicalPointType& metric, const ScalarType tolerance,
    const bool& check_orientation = true) {
  // -- Auxiliary data --
  constexpr std::size_t kParametricDimensions_ = parametric_dimension;
  constexpr auto subelement_vertex_ids =
      HyperCube<kParametricDimensions_>::SubElementVerticesToFace();
  constexpr auto opposite_face_list =
      HyperCube<kParametricDimensions_>::GetOppositeFaces();
  constexpr std::size_t number_of_element_faces = opposite_face_list.size();
  constexpr std::size_t number_of_element_vertices{algorithms::IntPower(
      static_cast<std::size_t>(2), kParametricDimensions_)};
  const std::size_t number_of_patches =
      corner_vertices.size() / number_of_element_vertices;

  // Consistency check
  if (!(corner_vertices.size() % number_of_element_vertices == 0)) {
    Logger::TerminatingError(
        "Number of corner vertices invalid. Must be a multiple of the number "
        "of vertices per patch");
  }

  // -- Determine connectivity --
  // Calculate face centers
  std::vector<PhysicalPointType> face_center_vertices;
  face_center_vertices.reserve(number_of_patches * number_of_element_faces);

  // Start the actual loop
  // (Instead of using the mean of the face vertices, using the sum)
  for (std::size_t i_spline{}; i_spline < number_of_patches; i_spline++) {
    const std::size_t patch_vertex_index_offset =
        i_spline * number_of_element_vertices;
    const std::size_t patch_face_index_offset =
        i_spline * number_of_element_faces;
    for (std::size_t i_face{}; i_face < number_of_element_faces; i_face++) {
      face_center_vertices.push_back(
          corner_vertices[patch_vertex_index_offset +
                          subelement_vertex_ids[i_face][0]]);
      for (std::size_t i_point{1}; i_point < subelement_vertex_ids[0].size();
           i_point++) {
        face_center_vertices[patch_face_index_offset + i_face] +=
            corner_vertices[patch_vertex_index_offset +
                            subelement_vertex_ids[i_face][i_point]];
      }
    }
  }

  // Return connectivity
  // Conditional is required for individual use without orientation for higher
  // dimensional or embedded types
  if (check_orientation) {
    return FindConnectivityFromCenters<kParametricDimensions_, true>(
        face_center_vertices, metric, tolerance);
  } else {
    return FindConnectivityFromCenters<kParametricDimensions_, false>(
        face_center_vertices, metric, tolerance);
  }
}

/**
 * @brief Finds duplicate Points and returns a list with indices that can be
 * used to build this list.
 *
 * example:
 * the list
 * [[0,1],[2,1],[0,1],[0,2]]
 * could return
 * [2,1,2,0] (order depends on orientation metric)
 */
template <std::size_t physical_dimension, typename ScalarType>
std::vector<std::size_t> IndexUniquePointList(
    const std::vector<Point<physical_dimension, ScalarType>>&
        original_point_list,
    const Point<physical_dimension, ScalarType> metric,
    const ScalarType tolerance = 1e-5) {
  Logger::Logging("Indexing unique point list");
  // Assure Metric is normed and non-zero
  if (metric.EuclidianNorm() <= 0) {
    Logger::TerminatingError("Metric is not normed or zero");
  }
  const Point<physical_dimension, ScalarType> normed_metric =
      metric * (static_cast<ScalarType>(1.) / metric.EuclidianNorm());

  // Store information in Auxiliary Values
  const std::size_t number_of_center_vertices{original_point_list.size()};
  const ScalarType tolerance_squared{tolerance * tolerance};

  // Init unique_indices and metric value
  // (-1 : untouched)
  // {in c++20 this expression could be constexpr}
  std::vector<std::size_t> unique_indices(
      // Size of vector
      number_of_center_vertices,
      // default value
      static_cast<std::size_t>(-1));

  // Initialize Metric
  std::vector<ScalarType> scalar_metric{};
  scalar_metric.reserve(number_of_center_vertices);

  // Check Metric Dimension and Vector Size
  for (unsigned int i{}; i < number_of_center_vertices; i++) {
    scalar_metric.push_back(normed_metric * original_point_list[i]);
  }

  // Sort Metric Vector
  const auto metric_order_indices = algorithms::IndexListSort(scalar_metric);

  // Start Uniquifying
  Logger::ExtendedInformation("Start unique indexing of control points");
  std::size_t number_of_new_points{};
  for (std::size_t lower_limit{0}; lower_limit < number_of_center_vertices - 1;
       lower_limit++) {
    // Point already processed
    if (unique_indices[metric_order_indices[lower_limit]] !=
        static_cast<std::size_t>(-1)) {
      continue;
    }

    // Point has not been processed -> add it to new point list
    unique_indices[metric_order_indices[lower_limit]] = number_of_new_points;

    // Now check allowed range for duplicates
    unsigned int upper_limit = lower_limit + 1;
    while (upper_limit < number_of_center_vertices &&
           scalar_metric[metric_order_indices[upper_limit]] -
                   scalar_metric[metric_order_indices[lower_limit]] <
               tolerance) {
      const bool found_duplicate =
          (original_point_list[metric_order_indices[lower_limit]] -
           original_point_list[metric_order_indices[upper_limit]])
              .SquaredEuclidianNorm() < tolerance_squared;
      if (found_duplicate) {
        if (unique_indices[metric_order_indices[upper_limit]] !=
            static_cast<std::size_t>(-1)) {
          Logger::TerminatingError(
              "Failure in indexing Unique Point List. "
              "Found more than two different indices in less than two "
              "tolerances proximity");
        }
        unique_indices[metric_order_indices[upper_limit]] =
            number_of_new_points;
      }
      upper_limit++;
    }
    number_of_new_points++;
  }

  // Special case
  const auto& last_index = metric_order_indices.size() - 1;
  if (unique_indices[metric_order_indices[last_index]] ==
      static_cast<std::size_t>(-1)) {
    unique_indices[metric_order_indices[last_index]] = number_of_new_points;
  }
  Logger::Logging("Found " + std::to_string(number_of_new_points) +
                  " unique points out of " +
                  std::to_string(number_of_center_vertices) + " points");
  return unique_indices;
}

/**!
 * Provide all necessary data for mfem export
 *
 * @attention This function is to be used with caution, as it makes several
 * assumptions about the underlying structure. It always assumes that the
 * parametric dimension matches the physical dimension of the problem
 *
 * @tparam PhysicalPointType array-type for coordinates
 * @tparam ScalarType
 * @param corner_vertices vector<PhysicalPointType> vector of corner vertices
 * @param metric used for preordering the vertices along a line
 * @param tolerance tolerance (distance between two vertices that are joined)
 * @return std::make_tuple(connectivity, vertex_ids, edge_information,
 *                        boundaries)
 */
template <typename PhysicalPointType,
          typename ScalarType = typename PhysicalPointType::ScalarType_>
auto ExtractMFEMInformation(
    const std::vector<PhysicalPointType>& corner_vertices,
    const PhysicalPointType& metric, const ScalarType& tolerance = 1e-5) {
  // Inform user in debug build with extensive export, that the function is
  // using strong assumptions
  Logger::Logging(
      "MFEM-data-extraction assumes Parametric dimension and spatial dimension "
      "to be equal");
  Logger::Logging("Start MFEM-data extraction");

  // -- Auxiliary data --
  constexpr std::size_t kParametricDimensions_ = PhysicalPointType{}.size();
  constexpr auto subelement_vertex_ids =
      HyperCube<kParametricDimensions_>::SubElementVerticesToFace();
  constexpr std::size_t number_of_element_vertices{algorithms::IntPower(
      static_cast<std::size_t>(2), kParametricDimensions_)};
  constexpr std::size_t number_of_vertices_per_face{algorithms::IntPower(
      static_cast<std::size_t>(2), kParametricDimensions_ - 1)};
  const std::size_t number_of_patches =
      corner_vertices.size() / number_of_element_vertices;
  constexpr auto faces_orthogonal_to_parametric_dimension =
      algorithms::HyperCube<
          kParametricDimensions_>::GetNormalFaceIndicesToParametricDimension();

  constexpr std::size_t number_of_element_faces_with_same_knotvectors{
      2 * kParametricDimensions_ - 2};

  // Consistency check
  if (corner_vertices.size() % number_of_element_vertices != 0) {
    Logger::TerminatingError(
        "Number of corner-vertices invalid. Maybe the spline does not fulfill "
        "paradim=dim");
  }

  // Find connectivity using face centers
  const auto connectivity = FindConnectivityFromCorners<kParametricDimensions_>(
      corner_vertices, metric, tolerance);
  Logger::ExtendedInformation("Connectivity determined");

  // -- Enumerate Vertices --
  const auto vertex_ids =
      algorithms::IndexUniquePointList(corner_vertices, metric);
  Logger::ExtendedInformation("Enumerated Vertices using uniquified points");

  // -- Determine and enumerate knot-vectors --// Enumerate knot vectors
  std::vector<std::array<std::size_t, kParametricDimensions_>> knot_vector_ids(
      number_of_patches,
      // Fill with array created by lambda function
      []() {
        std::array<std::size_t, kParametricDimensions_> init_{};
        init_.fill(static_cast<std::size_t>(-1));
        return init_;
      }());

  // Initialize knot_vector_counter
  std::size_t n_assigned_knot_vectors{};

  // index of splines which has to be control if they have neighbors
  std::vector<std::size_t> splines_to_go{};

  // Count Boundaries at the same time to preallocate boundary vector
  std::size_t number_of_boundaries{};

  // loop over all splines
  for (std::size_t i_spline{}; i_spline < number_of_patches; i_spline++) {
    // loop over all dimensions
    for (std::size_t i_dim{}; i_dim < kParametricDimensions_; i_dim++) {
      // check if spline has already a knot vector id
      if (knot_vector_ids[i_spline][i_dim] != static_cast<size_t>(-1)) {
        continue;
      }

      // insert index of start spline
      splines_to_go.push_back(i_spline);

      // retrieve relevant faces for propagation
      const auto& faces = faces_orthogonal_to_parametric_dimension[i_dim];

      // Assign knot vector id to current spline
      knot_vector_ids[i_spline][i_dim] = n_assigned_knot_vectors;

      // run through all splines until the vector splines_to_go is empty
      while (!splines_to_go.empty()) {
        // get last spline_id of the vector
        const std::size_t current_spline_id = splines_to_go.back();
        // delete the set current_spline_id from the vector
        splines_to_go.pop_back();

        // run through all neighbors
        for (std::size_t i_neighbor{};
             i_neighbor < number_of_element_faces_with_same_knotvectors;
             i_neighbor++) {
          const std::size_t& face_to_check = faces[i_neighbor];
          // Check if boundary
          if (connectivity[current_spline_id][face_to_check] ==
              static_cast<std::size_t>(-1)) {
            number_of_boundaries++;
            continue;
          }
          const std::size_t& current_neighbor =
              connectivity[current_spline_id][face_to_check];

          // if neighbor spline has the same  knot vector id as the
          // current spline continue --> Loop is closed
          if (knot_vector_ids[current_neighbor][i_dim] ==
              n_assigned_knot_vectors) {
            continue;
          }

          // Check if neighbor is different knot_vector_id, if true something
          // went wrong --> error
          if (knot_vector_ids[current_neighbor][i_dim] !=
              static_cast<std::size_t>(-1)) {
            Logger::TerminatingError(
                "Connectivity error. Knotvector mismatch between two "
                "neighboring splines");
          }

          // set new knot vector
          knot_vector_ids[current_neighbor][i_dim] = n_assigned_knot_vectors;
          // Add neighbor spline to splines_to_go
          splines_to_go.push_back(current_neighbor);
        }
      }
      n_assigned_knot_vectors++;
    }
  }

  // The boundary count is only valid for 2D elements, as all faces correspond
  // to exactly one knot_vector dimension. In higher order problems, like 3D
  // the same face needs to be checked for n_dim-1 knot_vector dimensions, the
  // number of boundaries therefore needs to be updated accordingly
  assert(number_of_boundaries % (kParametricDimensions_ - 1) == 0);
  number_of_boundaries /= (kParametricDimensions_ - 1);

  // Lastly the knotvectors need exporting only once per edge

  // map for the edges -> key is the edge with the vertices | value is the
  // knot vector id
  std::map<std::pair<std::size_t, std::size_t>, std::size_t>
      edge_uniquify_map{};

  constexpr auto local_face_vertex_index =
      algorithms::HyperCube<kParametricDimensions_>::EdgeVertexIndices();
  constexpr std::size_t number_of_edges_per_element =
      local_face_vertex_index.size();

  for (std::size_t i_spline{}; i_spline < number_of_patches; i_spline++) {
    for (std::size_t i_edge{}; i_edge < number_of_edges_per_element; i_edge++) {
      // Edges are ordered by dimension using integer division
      std::size_t i_dim{i_edge / (kParametricDimensions_ == 2 ? 2 : 4)};
      // Add element to edge
      edge_uniquify_map[
          // Create Pair of vertex ids as key
          std::make_pair(
              // First Vertex
              vertex_ids[i_spline * number_of_element_vertices +
                         local_face_vertex_index[i_edge][0]],
              // Second Vertex
              vertex_ids[i_spline * number_of_element_vertices +
                         local_face_vertex_index[i_edge][1]])] =
          // set the value
          knot_vector_ids[i_spline][i_dim];
    }
  }

  // Transform data into vector
  std::vector<std::array<std::size_t, 3>> edge_information;
  edge_information.reserve(edge_uniquify_map.size());

  for (const auto& [edge_vertices, i_knot_id] : edge_uniquify_map) {
    edge_information.push_back(
        {i_knot_id, edge_vertices.first, edge_vertices.second});
  }
  Logger::ExtendedInformation("Knotvector and Edge assignment successful");

  // Last step is to identify all boundary elements
  std::vector<std::array<std::size_t, number_of_vertices_per_face>> boundaries;
  boundaries.reserve(number_of_boundaries);
  // @todo Provide functions to set appropriate boundary id here
  for (std::size_t i_spline{}; i_spline < number_of_patches; i_spline++) {
    for (std::size_t i_face{}; i_face < 2 * kParametricDimensions_; i_face++) {
      if (connectivity[i_spline][i_face] == static_cast<std::size_t>(-1)) {
        // write the vertices
        const auto& local_face_vertices = subelement_vertex_ids[i_face];
        std::array<std::size_t, number_of_vertices_per_face> boundary_element{};
        for (std::size_t i_local_vertex_on_face{};
             i_local_vertex_on_face < local_face_vertices.size();
             i_local_vertex_on_face++) {
          boundary_element[i_local_vertex_on_face] =
              vertex_ids[i_spline * number_of_element_vertices +
                         local_face_vertices[i_local_vertex_on_face]];
        }
        boundaries.push_back(boundary_element);
      }
    }
  }

  // Put all data into a tuple for easy access
  return std::make_tuple(connectivity, vertex_ids, edge_information,
                         boundaries);
}

}  // namespace bezman::utils::algorithms
#endif  // UTILS_UNIQUIFY_POINT_UNIQUIFIER_HPP
