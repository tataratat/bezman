#ifndef SRC_UTILS_EXPORT_INC
#define SRC_UTILS_EXPORT_INC

std::string Export::ensureFileExtension(const std::string& filename,
                                        const std::string& obj_extension) {
  std::size_t ext_start = filename.find_last_of(".");
  const std::string extension = filename.substr(ext_start + 1);
  if (extension == obj_extension) {
    return filename;
  } else {
    return (filename + ".") + obj_extension;
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2IRITfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file) {
  // Irit only support 9 dimensional spaces
  static_assert((PhysicalPointType::kSpatialDimension < 10),
                "IRIT natively only supports up to 9D spaces");

  // Increase user readability
  const std::string indentation{"\t"};

  // Write Header
  export_file << indentation << indentation << "[MULTIVAR BEZIER "
              << parametric_dimension << " ";

  // Write different degrees
  for (std::size_t i{}; i < parametric_dimension; i++) {
    // IRIT expexts order not degrees
    export_file << spline.GetDegrees()[i] + 1 << " ";
  }

  // Point Type
  export_file << "E" << PhysicalPointType::kSpatialDimension << "\n";

  // Write All Control Points
  for (std::size_t i_point{0}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    export_file << indentation << indentation << indentation
                << spline.control_points[i_point] << "\n";
  }

  // Write Footer
  export_file << indentation << indentation << "]";
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2JSONfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file, const bool base64encoding) {
  // Increase user readability
  const std::string indentation{"      "};

  // Write Dimensionality
  export_file << indentation
              << "\"dim\" : " << PhysicalPointType::kSpatialDimension << ",\n";

  // Write Parametric Dimensionality
  export_file << indentation << "\"para_dim\" : " << parametric_dimension
              << ",\n";

  // Write degrees
  export_file << indentation << "\"degrees\" : [" << spline.GetDegrees()[0];
  for (std::size_t i{1}; i < parametric_dimension; i++) {
    // IRIT expexts order not degrees
    export_file << ", " << spline.GetDegrees()[i];
  }
  export_file << "],\n";

  // Write All Control Points
  if (base64encoding) {
    export_file << indentation << "\"control_points\" : \""
                << utils::Base64::Encode(spline.control_points) << "\"\n";
  } else {
    export_file << indentation << "\"control_points\" : [\n"
                << indentation << "  " << spline.control_points[0];
    for (std::size_t i_point{1}; i_point < spline.NumberOfControlPoints;
         i_point++) {
      export_file << ",\n"
                  << indentation << "  " << spline.control_points[i_point];
    }
    export_file << "\n" << indentation << "]\n";
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::format2XMLfile(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    std::ofstream& export_file) {
  // Handy alias
  using BezierSplineType =
      BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;
  // Increase user readability
  const std::string indentation{"  "};

  // Write Header
  export_file << indentation << "<SplineEntry splDim=\"" << parametric_dimension
              << "\" spaceDim=\"" << PhysicalPointType::kSpatialDimension
              << "\" numOfCntrlPntVars=\""
              << PhysicalPointType ::kSpatialDimension << "\" numCntrlPnts=\""
              << spline.NumberOfControlPoints << "\" numOfEleVars=\"" << 0
              << "\">\n";

  // Write ControlPoint Variable Names
  export_file << indentation << "<cntrlPntVarNames>\n";
  for (std::size_t i{}; i < PhysicalPointType ::kSpatialDimension; i++) {
    export_file << indentation << "  x" << i << " ";
  }
  export_file << "\n" << indentation << "</cntrlPntVarNames>\n";

  // Write Points
  export_file << indentation << "<cntrlPntVars>";
  for (std::size_t i_point{}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    // Spacing
    export_file << "\n" << indentation << indentation;
    for (std::size_t i_dim{}; i_dim < PhysicalPointType ::kSpatialDimension;
         i_dim++) {
      export_file << std::setw(12) << std::setprecision(10)
                  << spline.control_points[i_point][i_dim] << indentation;
    }
  }
  export_file << "\n" << indentation << "</cntrlPntVars>\n";

  // Write Weights
  export_file << indentation << "<wght>\n" << indentation << indentation;
  for (std::size_t i_point{}; i_point < spline.NumberOfControlPoints;
       i_point++) {
    export_file << "1.0  ";  // Add \n if required
  }
  export_file << "\n" << indentation << "</wght>\n";

  // Write Degrees
  export_file << indentation << "<deg>\n" << indentation << indentation;
  for (std::size_t i_par_dim{};
       i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
    export_file << spline.GetDegrees()[i_par_dim]
                << indentation;  // Add \n if required
  }
  export_file << "\n" << indentation << "</deg>\n";

  // Knotvectors (required)
  export_file << indentation << "<kntVecs>\n";
  for (std::size_t i_par_dim{};
       i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
    export_file << indentation << indentation << "<kntVec>\n"
                << indentation << indentation;
    // Write All Control Points
    for (std::size_t i_point{0};
         i_point < 2 * spline.GetDegrees()[i_par_dim] + 2; i_point++) {
      export_file << indentation
                  << (i_point < (spline.GetDegrees()[i_par_dim] + 1) ? "0."
                                                                     : "1.");
    }
    export_file << "\n" << indentation << indentation << "</kntVec>\n";
  }
  export_file << indentation << "</kntVecs>\n";

  // Write Footer
  export_file << indentation << "</SplineEntry>\n";
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::GuessByExtension(
    const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>&
        spline,
    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  GuessByExtension(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::GuessByExtension(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // In case of a absolute path strip filename
  std::size_t last_found_id = filename.find_last_of("/\\");
  const std::string file_without_path = filename.substr(last_found_id + 1);

  // Check file extension
  std::size_t extension_start_id = file_without_path.find_last_of(".");
  const std::string extension_str =
      file_without_path.substr(extension_start_id + 1);

  if (extension_str == "itd") {
    // Export as IRIT
    AsIRIT(spline_group, filename);
  } else if (extension_str == "xml") {
    // Export as XML
    AsXML(spline_group, filename);
  } else if (extension_str == "json") {
    // Export as custom JSON
    AsJSON(spline_group, filename);
  } else if (extension_str == "mesh") {
    // Export as MFEM
    if constexpr (parametric_dimension == 2 &&
                  PhysicalPointType::kSpatialDimension == 2) {
      AsMFEM(spline_group, filename);
    } else {
      assert(false);
    }
  } else {
    throw std::invalid_argument("File Extension not supported");
  }
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsIRIT(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsIRIT(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsIRIT(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "itd"));

  // Write Header
  export_file << "[OBJECT BEZIERGROUP";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // Write Header
    export_file << "\n\t[OBJECT BEZIERSPLINE" << i;
    // Small spacer for Readability
    export_file << "\n";
    // Export individual splines seperately
    format2IRITfile(spline_group[i], export_file);
    // Finalize Current Spline
    export_file << "\n\t]";
  }

  // Write Footer
  export_file << "\n]";
  export_file.close();
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsXML(const BezierSpline<parametric_dimension, PhysicalPointType,
                                      ScalarType>& spline,
                   const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsXML(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsXML(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "xml"));

  // Write Header
  export_file << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n";

  // Start Spline List
  export_file << "<SplineList SplineType = \"1\" NumberOfSplines =\""
              << spline_group.size() << "\">\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // Export individual splines seperately
    format2XMLfile(spline_group[i], export_file);
  }

  // Write Footer
  export_file << "</SplineList>\n";
  export_file.close();
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename, const bool base64encoding) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsJSON(spline_group, filename, base64encoding);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsJSON(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename, const bool base64encoding) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "json"));

  // Write Header
  export_file << "{\n  \"Name\" : \"BezierSplineGroup\",\n";

  // Number of Splines
  export_file << "  \"NumberOfSplines\" : " << spline_group.size() << ",\n";

  // Encoding
  export_file << "  \"Base64Encoding\" : " << base64encoding << ",\n";

  // Start Spline Export
  export_file << "  \"SplineList\" : [\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // ID and Type are set here
    export_file
        << "    {\n      \"SplineType\" : \"Bezier\",\n      \"SplineID\" : "
        << i << ",\n";
    // Export individual splines seperately
    format2JSONfile(spline_group[i], export_file, base64encoding);
    export_file << "    }";
    if (i < spline_group.size() - 1) {
      export_file << ",\n";
    } else {
      export_file << "\n";
    }
  }

  // Write Footer
  export_file << "  ]\n}";
  export_file.close();
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsMFEM(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType>& spline,
                    const std::string& filename) {
  BezierSplineGroup<parametric_dimension, PhysicalPointType, ScalarType>
      spline_group{spline};
  AsMFEM(spline_group, filename);
}

template <std::size_t parametric_dimension, typename PhysicalPointType,
          typename ScalarType>
void Export::AsMFEM(
    const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                            ScalarType>& spline_group,
    const std::string& filename) {
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "mesh"));

  // MFEM format currently restricted to 2D2D and 3D3D
  static_assert(parametric_dimension == PhysicalPointType::kSpatialDimension,
                "Mismatch between parametric and physical dimension");

  // Retrieve Connectivity
  const auto connectivity =
      uniquify::GetConnectivityForSplineGroup(spline_group);
  // @todo throw exception
  assert(connectivity.size() > 0);

  /*
   * Further steps:
   * 1. Calculate Vertices based on connectivity.
   *    - Create a new matrix that contains vertex-ids
   *    - Loop over elements and edges (increasing ekement id)
   *    - Add Vertex ID in two steps
   *      1. Retrieve all neighbors, that have a lower element-ID
   *      2. Use their existing Vertex IDs on connected faces
   *      3. Assign new Ids to unset vertices
   *      - As the element ID is increasing, no double association will take
   *        place
   * 2. Determine KnotVector IDs
   *    - Create a new Matrix, that has the same dimensions as the connectivity
   *    - Assign all Entries in new matrix to a dummy-value, that represents
   *      "untouched", e.g., -1
   *    - Pseudoalgorithm:
   *      for i_el in elements:
   *        for i_face in face:
   *          if knotV[i_el, i_face] is touched:
   *            continue
   *          else:
   *            -Assign new knotV ID to this face
   *            while True:
   *              - Get next element that connects to this face
   *              - Assign new knotVID
   *              if connected face is (a) already assigned this id or (b)
   *                            boundary:
   *                break
   *              else if connected face is already assigned a different id
   *                -> ERROR
   *            - Do the same while loop in other direction
   * 3. Export data into mesh format
   */

  // Auxiliary values
  const std::size_t number_of_splines{spline_group.size()};
  constexpr auto face_vertex_indices =
      algorithms::HyperCube<parametric_dimension>::FaceVertexIndices();
  constexpr auto opposite_faces =
      algorithms::HyperCube<parametric_dimension>::GetOppositeFaces();
  constexpr std::size_t number_of_element_faces{opposite_faces.size()};
  constexpr std::size_t number_of_vertices_per_face{
      face_vertex_indices[0].size()};
  constexpr std::size_t number_of_element_vertices{
      algorithms::IntPower(static_cast<std::size_t>(2), parametric_dimension)};

  int number_of_boundaries = 0;

  // Start with 1. (overflow as -1 to imply untouched)
  std::vector<std::array<std::size_t, number_of_element_vertices>> vertex_id(
      spline_group.size(),
      // Fill with array created by lambda function
      []() {
        std::array<std::size_t, number_of_element_vertices> a{};
        a.fill(static_cast<std::size_t>(-1));
        return a;
      }());

  // Start loop
  std::size_t n_assigned_vertices{};
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    for (std::size_t i_face{}; i_face < number_of_element_faces; i_face++) {
      // Check if neighbor id is smaller than current spline id
      if (connectivity[i_spline][i_face] < i_spline) {
        const std::size_t& face_in_neighbor_spline = opposite_faces[i_face];
        const std::size_t& neighbor_id = connectivity[i_spline][i_face];
        // This is a bit inefficient, because some vertices might be touched
        // multiple times
        for (std::size_t i_vertex_in_face{};
             i_vertex_in_face < number_of_vertices_per_face;
             i_vertex_in_face++) {
          vertex_id[i_spline][face_vertex_indices[i_face][i_vertex_in_face]] =
              vertex_id[neighbor_id]
                       [face_vertex_indices[face_in_neighbor_spline]
                                           [i_vertex_in_face]];
        }
      }
    }
    // Enumerate remaining vertices
    for (std::size_t i_vertex{}; i_vertex < number_of_element_vertices;
         i_vertex++) {
      if (vertex_id[i_spline][i_vertex] == static_cast<std::size_t>(-1)) {
        vertex_id[i_spline][i_vertex] = n_assigned_vertices;
        n_assigned_vertices++;
      }
    }
  }

  // Continue with 2.
  // @todo Replace all the assertions with exceptions
  // Enumerate knot vectors
  std::vector<std::array<std::size_t, number_of_element_faces>> knot_vector_ids(
      number_of_splines,
      // Fill with array created by lambda function
      []() {
        std::array<std::size_t, number_of_element_faces> init_{};
        init_.fill(static_cast<std::size_t>(-1));
        return init_;
      }());

  // Start the loop over all knot_vectors
  std::size_t n_assigned_knot_vector{};
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    for (std::size_t i_face{}; i_face < number_of_element_faces; i_face++) {
      // Skip knot vector if assigned
      if (knot_vector_ids[i_spline][i_face] != static_cast<std::size_t>(-1)) {
        continue;
      } else {
        knot_vector_ids[i_spline][i_face] = n_assigned_knot_vector;
        // Initialize Variables
        std::size_t neighbor_element = i_spline;
        std::size_t neighbor_face = i_face;
        // Loop into the direction of the neighboring spline
        while (connectivity[neighbor_element][neighbor_face] != -1) {
          neighbor_element = connectivity[i_spline][i_face];
          neighbor_face = opposite_faces[i_face];
          // If the loop is closed -> stop
          if (knot_vector_ids[neighbor_element][neighbor_face] ==
              n_assigned_knot_vector) {
            break;
          }
          // Check if there might be an orientation error
          assert(knot_vector_ids[neighbor_element][neighbor_face] ==
                 static_cast<std::size_t>(-1));
          // Assign knot vector
          knot_vector_ids[neighbor_element][neighbor_face] =
              n_assigned_knot_vector;
          // Update neighbor face
          neighbor_face = opposite_faces[neighbor_face];
          // If the loop is closed -> stop
          if (knot_vector_ids[neighbor_element][neighbor_face] ==
              n_assigned_knot_vector) {
            break;
          }
          assert(knot_vector_ids[neighbor_element][neighbor_face] ==
                 static_cast<std::size_t>(-1));
          knot_vector_ids[neighbor_element][neighbor_face] =
              n_assigned_knot_vector;
        }
        // Continue the same loop into the other direction
        neighbor_face = opposite_faces[i_face];
        neighbor_element = i_spline;
        // Check if there might be an orientation error
        assert((knot_vector_ids[neighbor_element][neighbor_face] ==
                static_cast<std::size_t>(-1)) ||
               (knot_vector_ids[neighbor_element][neighbor_face] ==
                n_assigned_knot_vector));
        // Assign and loop in the other direction
        knot_vector_ids[neighbor_element][neighbor_face] =
            n_assigned_knot_vector;
        while (connectivity[neighbor_element][neighbor_face] != -1) {
          // Continue to next element if applies
          neighbor_element = connectivity[neighbor_element][neighbor_face];
          neighbor_face = opposite_faces[neighbor_face];
          // If the loop is closed -> stop
          if (knot_vector_ids[neighbor_element][neighbor_face] ==
              n_assigned_knot_vector) {
            break;
          }
          // Check if there might be an orientation error
          assert(knot_vector_ids[neighbor_element][neighbor_face] ==
                 static_cast<std::size_t>(-1));
          // Assign knot vector
          knot_vector_ids[neighbor_element][neighbor_face] =
              n_assigned_knot_vector;
          // Update neighbor face
          neighbor_face = opposite_faces[neighbor_face];
          // If the loop is closed -> stop
          if (knot_vector_ids[neighbor_element][neighbor_face] ==
              n_assigned_knot_vector) {
            break;
          }
          assert(knot_vector_ids[neighbor_element][neighbor_face] ==
                 static_cast<std::size_t>(-1));
          knot_vector_ids[neighbor_element][neighbor_face] =
              n_assigned_knot_vector;
        }
        // Update counter
        n_assigned_knot_vector++;
      }
    }
  }
  // Now all information is available... Write data into file accoring to
  // @deltaluki

  // Write Header
  export_file << "MFEM NURBS mesh v1.0\n";

  // Write dimension
  export_file << "\ndimension \n" << parametric_dimension << std::endl;

  // Write segment element
  export_file << "\nelements\n";
  export_file << number_of_splines << std::endl;

    for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
        // write attribute and geometry type TODO: use var for geo type and attr
        export_file << 1 << " " << 3;
        for (std::size_t i_vertex{};
             i_vertex < number_of_element_vertices; i_vertex++) {
            export_file << " " <<vertex_id[i_spline][i_vertex];
        }
        export_file << std::endl;
    }

    //TODO maybe a better solution
    // is needed in the vertices segment
    int max_vertices = 0;

    //write boundary segment
    int counter = 0;
    export_file << "\n\nboundary\n";
    export_file << "100\n"; //TODO working for solution
    for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
        auto element = connectivity[i_spline];

        for(int i = 0; i < 4; i++){
            //search max vertices
            if(element[i] > max_vertices){
                max_vertices = element[i];
            }
            if(element[i] == -1){
                // 1 and 3 edges (__)
                if(i == 0 || i == 3){
                    export_file << 1 << " " << 1 << " " << vertex_id[i_spline][i];
                    if(i == 4){
                        export_file << " "<< vertex_id[i_spline][0] << std::endl;
                    }else{
                        export_file << " " << vertex_id[i_spline][i+1] << std::endl;
                    }
                // 2 and 4 edges (||)
                }else if(i == 1 || i == 2){
                    export_file << 2 << " " << 1 << " " <<
                    vertex_id[i_spline][i] << " " <<vertex_id[i_spline][i+1]<< std::endl;
                }
            }
        }
    }

    // write edge boundary
    export_file << "\n\nedges\n";
    export_file << (number_of_splines * 4) << std::endl;
    for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
        for(int j = 0; j < 4; j++){
            int index_vertex_2 = j + 1;
            if(j == 3){
                index_vertex_2 = 0;
            }
            export_file << knot_vector_ids[i_spline][j] <<
            " " << vertex_id[i_spline][j] << " " << vertex_id[i_spline][index_vertex_2] << std::endl;
        }
    }
    export_file << std::endl;

    // write vertices
    export_file << "\n\nvertices\n";
    export_file << max_vertices + 1 << std::endl;

    using BezierSplineType =
            BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>;

    //write patches
    export_file << "\n\npatches\n";
    for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {

        // write knotvectors
        export_file << "\nknotvectors" << std::endl;
        export_file << BezierSplineType::kParametricDimensions << std::endl;

        //copy code from format2XMLfile  ;)
        for (std::size_t i_par_dim{};
             i_par_dim < BezierSplineType::kParametricDimensions; i_par_dim++) {
            // Write All Control Points
            export_file << spline_group[i_spline].GetDegrees()[i_par_dim] << " "
            << spline_group[i_spline].GetDegrees()[i_par_dim]+1 << " ";
            for (std::size_t i_point{0};
                 i_point < 2 * spline_group[i_spline].GetDegrees()[i_par_dim] + 2; i_point++) {
                export_file <<
                (i_point < (spline_group[i_spline].GetDegrees()[i_par_dim] + 1) ? "0.0 "
                                                                               : "1.0 ");
            }
            export_file << std::endl;
        }

        // write dimension
        export_file << "\ndimension" << std::endl;
        //TODO is this the correct value?
        export_file << BezierSplineType::kParametricDimensions << std::endl;

        //controlpoints_cartesian
        export_file << "\n\ncontrolpoints_cartesian\n";
        for (std::size_t i_point{}; i_point < spline_group[i_spline].NumberOfControlPoints;
             i_point++) {
            for (std::size_t i_dim{}; i_dim < PhysicalPointType ::kSpatialDimension;
                 i_dim++) {
                export_file << spline_group[i_spline].control_points[i_point][i_dim]
                            << " ";
            }
            //write weight
            export_file << " " << 1 << std::endl;
        }

        export_file << std::endl;

    }

  // Write Footer
  export_file.close();
}

#endif  // SRC_UTILS_EXPORT_INC
