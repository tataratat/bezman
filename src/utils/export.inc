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

template <typename SplineType>
void Export::format2IRITfile(const SplineType& spline,
                             std::ofstream& export_file) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  // Irit only support 9 dimensional spaces
  static_assert((SplineType::PhysicalPointType_::kSpatialDimension < 10),
                "IRIT natively only supports up to 9D spaces");

  // Increase user readability
  const std::string indentation{"\t"};

  // Write Header
  export_file << indentation << indentation << "[MULTIVAR BEZIER "
              << SplineType::kParametricDimensions << " ";

  // Write different degrees
  for (std::size_t i{}; i < SplineType::kParametricDimensions; i++) {
    // IRIT expects order not degrees
    export_file << spline.GetDegrees()[i] + 1 << " ";
  }

  // Point Type
  export_file << (isRationalBezierSpline_v<SplineType> ? "P" : "E")
              << SplineType::PhysicalPointType_::kSpatialDimension << "\n";

  // Write All Control Points
  for (std::size_t i_point{0}; i_point < spline.GetNumberOfControlPoints();
       i_point++) {
    export_file << indentation << indentation << indentation;
    if constexpr (!isRationalBezierSpline_v<SplineType>) {
      export_file << spline.control_points[i_point] << "\n";
    } else {
      export_file << "[" << spline.GetWeights()[i_point];
      for (std::size_t i_dim{0};
           i_dim < SplineType::PhysicalPointType_::kSpatialDimension; i_dim++) {
        export_file << " " << spline.GetWeightedControlPoints()[i_point][i_dim];
      }
      export_file << "]\n";
    }
  }

  // Write Footer
  export_file << indentation << indentation << "]";
}

template <typename SplineType>
void Export::format2JSONfile(const SplineType& spline,
                             std::ofstream& export_file,
                             const bool base64encoding) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  // Increase user readability
  const std::string indentation{"      "};

  // Write Dimensionality
  export_file << indentation << "\"dim\" : "
              << SplineType::PhysicalPointType_::kSpatialDimension << ",\n";

  // Write Parametric Dimensionality
  export_file << indentation
              << "\"para_dim\" : " << SplineType::kParametricDimensions
              << ",\n";

  // Write degrees
  export_file << indentation << "\"degrees\" : [" << spline.GetDegrees()[0];
  for (std::size_t i{1}; i < SplineType::kParametricDimensions; i++) {
    // IRIT expexts order not degrees
    export_file << ", " << spline.GetDegrees()[i];
  }
  export_file << "],\n";

  // Write All Control Points
  if constexpr (!isRationalBezierSpline_v<SplineType>) {
    if (base64encoding) {
      export_file << indentation << "\"control_points\" : \""
                  << utils::Base64::Encode(spline.control_points) << "\"\n";
    } else {
      export_file << indentation << "\"control_points\" : [\n"
                  << indentation << "  " << spline.control_points[0];
      for (std::size_t i_point{1}; i_point < spline.GetNumberOfControlPoints();
           i_point++) {
        export_file << ",\n"
                    << indentation << "  " << spline.control_points[i_point];
      }
      export_file << "\n" << indentation << "]\n";
    }
  } else {
    // JSON expects non-weighted control-points
    std::vector<typename SplineType::PhysicalPointType_> ctps{};
    ctps.reserve(spline.GetNumberOfControlPoints());
    for (std::size_t i_point{}; i_point < spline.GetNumberOfControlPoints();
         i_point++) {
      ctps.push_back(spline.GetWeightedControlPoints()[i_point] /
                     spline.GetWeights()[i_point]);
    }
    if (base64encoding) {
      export_file << indentation << "\"control_points\" : \""
                  << utils::Base64::Encode(ctps) << "\",\n";
      export_file << indentation << "\"weights\" : \""
                  << utils::Base64::Encode(spline.GetWeights()) << "\"\n";
    } else {
      // Write Control Points
      export_file << indentation << "\"control_points\" : [\n"
                  << indentation << "  " << ctps[0];
      for (std::size_t i_point{1}; i_point < spline.GetNumberOfControlPoints();
           i_point++) {
        export_file << ",\n" << indentation << "  " << ctps[i_point];
      }
      export_file << "\n" << indentation << "],\n";
      // Write Weights
      export_file << indentation << "\"weights\" : [\n"
                  << indentation << "  " << spline.GetWeights()[0];
      for (std::size_t i_point{1}; i_point < spline.GetNumberOfControlPoints();
           i_point++) {
        export_file << ",\n"
                    << indentation << "  " << spline.GetWeights()[i_point];
      }
      export_file << "\n" << indentation << "]\n";
    }
  }
}

template <typename SplineType>
void Export::format2XMLfile(const SplineType& spline,
                            std::ofstream& export_file) {
  // @TODO: Implement for rational types
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(isBezierSpline_v<SplineType>,
                "Implemented only for bezier spline-types.");

  // Increase user readability
  const std::string indentation{"  "};

  // Write Header
  export_file << indentation << "<SplineEntry splDim=\""
              << SplineType::kParametricDimensions << "\" spaceDim=\""
              << SplineType::PhysicalPointType_::kSpatialDimension
              << "\" numOfCntrlPntVars=\""
              << SplineType::PhysicalPointType_ ::kSpatialDimension
              << "\" numCntrlPnts=\"" << spline.GetNumberOfControlPoints()
              << "\" numOfEleVars=\"" << 0 << "\">\n";

  // Write ControlPoint Variable Names
  export_file << indentation << "<cntrlPntVarNames>\n";
  for (std::size_t i{}; i < SplineType::PhysicalPointType_ ::kSpatialDimension;
       i++) {
    export_file << indentation << "  x" << i << " ";
  }
  export_file << "\n" << indentation << "</cntrlPntVarNames>\n";

  // Write Points
  export_file << indentation << "<cntrlPntVars>";
  for (std::size_t i_point{}; i_point < spline.GetNumberOfControlPoints();
       i_point++) {
    // Spacing
    export_file << "\n" << indentation << indentation;
    for (std::size_t i_dim{};
         i_dim < SplineType::PhysicalPointType_ ::kSpatialDimension; i_dim++) {
      export_file << std::setw(12) << std::setprecision(10)
                  << spline.control_points[i_point][i_dim] << indentation;
    }
  }
  export_file << "\n" << indentation << "</cntrlPntVars>\n";

  // Write Weights
  export_file << indentation << "<wght>\n" << indentation << indentation;
  for (std::size_t i_point{}; i_point < spline.GetNumberOfControlPoints();
       i_point++) {
    export_file << "1.0  ";  // Add \n if required
  }
  export_file << "\n" << indentation << "</wght>\n";

  // Write Degrees
  export_file << indentation << "<deg>\n" << indentation << indentation;
  for (std::size_t i_par_dim{}; i_par_dim < SplineType::kParametricDimensions;
       i_par_dim++) {
    export_file << spline.GetDegrees()[i_par_dim]
                << indentation;  // Add \n if required
  }
  export_file << "\n" << indentation << "</deg>\n";

  // Knotvectors (required)
  export_file << indentation << "<kntVecs>\n";
  for (std::size_t i_par_dim{}; i_par_dim < SplineType::kParametricDimensions;
       i_par_dim++) {
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

template <typename SplineType>
void Export::GuessByExtension(const SplineType& spline,
                              const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  BezierGroup<SplineType> spline_group{spline};
  GuessByExtension(spline_group, filename);
}

template <typename SplineType>
void Export::GuessByExtension(const BezierGroup<SplineType>& spline_group,
                              const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  // In case of a absolute path strip filename
  std::size_t last_found_id = filename.find_last_of("/\\");
  const std::string file_without_path = filename.substr(last_found_id + 1);

  // Check file extension
  std::size_t extension_start_id = file_without_path.find_last_of(".");
  const std::string extension_str =
      file_without_path.substr(extension_start_id + 1);

  if (extension_str == "itd") {
    Logger::Logging("Identified IRIT export format");
    // Export as IRIT
    AsIRIT(spline_group, filename);
  } else if (extension_str == "xml") {
    Logger::Logging("Identified XML export format");
    // Export as XML
    AsXML(spline_group, filename);
  } else if (extension_str == "json") {
    Logger::Logging("Identified json export format");
    Logger::ExtendedInformation("Json Export format uses base64 as default");
    Logger::ExtendedInformation("Use AsJSON function to change to text output");
    // Export as custom JSON
    AsJSON(spline_group, filename);
  } else if (extension_str == "mesh") {
    Logger::Logging("Identified MFEM export format");
    // Export as MFEM
    if constexpr (SplineType::kParametricDimensions ==
                      SplineType::PhysicalPointType_::kSpatialDimension &&
                  (SplineType::kParametricDimensions ==
                       static_cast<std::size_t>(2) ||
                   SplineType::kParametricDimensions ==
                       static_cast<std::size_t>(3))) {
      AsMFEM(spline_group, filename);
    } else {
      Logger::Error("MFEM only supports 2D->2D or 3D->3D Splines");
      Logger::TerminatingError("Not Implemented");
    }
  } else {
    throw std::invalid_argument("File Extension not supported");
  }
}

template <typename SplineType>
void Export::AsIRIT(const SplineType& spline, const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  BezierGroup<SplineType> spline_group{spline};
  AsIRIT(spline_group, filename);
}

template <typename SplineType>
void Export::AsIRIT(const BezierGroup<SplineType>& spline_group,
                    const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  Logger::UserInfo("Starting IRIT Export as : " + filename);
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
    // Export individual splines separately
    format2IRITfile(spline_group[i], export_file);
    // Finalize Current Spline
    export_file << "\n\t]";
  }

  // Write Footer
  export_file << "\n]";
  export_file.close();
  Logger::UserInfo("Finished IRIT Export");
}

template <typename SplineType>
void Export::AsXML(const SplineType& spline, const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  BezierGroup<SplineType> spline_group{spline};
  AsXML(spline_group, filename);
}

template <typename SplineType>
void Export::AsXML(const BezierGroup<SplineType>& spline_group,
                   const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  Logger::UserInfo("Starting XML Export as : " + filename);
  Logger::Warning(
      "XML Export inefficient (export as NURBS) and not standardized.");
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
    // Export individual splines separately
    format2XMLfile(spline_group[i], export_file);
  }

  // Write Footer
  export_file << "</SplineList>\n";
  export_file.close();
  Logger::UserInfo("Finished XML Export");
}

template <typename SplineType>
void Export::AsJSON(const SplineType& spline, const std::string& filename,
                    const bool base64encoding) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  BezierGroup<SplineType> spline_group{spline};
  AsJSON(spline_group, filename, base64encoding);
}

template <typename SplineType>
void Export::AsJSON(const BezierGroup<SplineType>& spline_group,
                    const std::string& filename, const bool base64encoding) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  Logger::UserInfo("Starting JSON Export as : " + filename);
  Logger::Logging(base64encoding ? "Exporting using Base64 Encoding"
                                 : "Exporting as Text-Only");
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "json"));

  // Write Header
  export_file << "{\n  \"Name\" : \"BezierGroup\",\n";

  // Number of Splines
  export_file << "  \"NumberOfSplines\" : " << spline_group.size() << ",\n";

  // Encoding
  export_file << "  \"Base64Encoding\" : " << base64encoding << ",\n";

  // Start Spline Export
  export_file << "  \"SplineList\" : [\n";

  // Write all Splines from Group into new file
  for (std::size_t i{0}; i < spline_group.size(); i++) {
    // ID and Type are set here
    export_file << "    {\n      \"SplineType\" : \""
                << (isRationalBezierSpline_v<SplineType> ? "RationalBezier"
                                                         : "Bezier")
                << "\",\n      \"SplineID\" : " << i << ",\n";
    // Export individual splines separately
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
  Logger::UserInfo("Finished XML Export");
}

template <typename SplineType>
void Export::AsMFEM(const SplineType& spline, const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  BezierGroup<SplineType> spline_group{spline};
  AsMFEM(spline_group, filename);
}

template <typename SplineType>
void Export::AsMFEM(const BezierGroup<SplineType>& spline_group,
                    const std::string& filename) {
  // Ensure Call from BezierType
  using namespace utils::type_traits;
  static_assert(
      isBezierSpline_v<SplineType> || isRationalBezierSpline_v<SplineType>,
      "Implemented only for (rational) bezier spline-types.");

  Logger::UserInfo("Starting MFEM Export as : " + filename);
  Logger::Warning(
      "MFEM Exports of derivatives will result in different enumerations!");
  // Create File Stream
  std::ofstream export_file;
  export_file.open(ensureFileExtension(filename, "mesh"));

  // -- Auxiliary Values --
  const std::size_t number_of_splines{spline_group.size()};
  constexpr std::size_t number_of_element_vertices{algorithms::IntPower(
      static_cast<std::size_t>(2), SplineType::kParametricDimensions)};

  // MFEM format currently restricted to 2D2D and 3D3D
  static_assert(SplineType::kParametricDimensions ==
                    SplineType::PhysicalPointType_::kSpatialDimension,
                "Mismatch between parametric and physical dimension");
  static_assert(
      SplineType::kParametricDimensions == static_cast<std::size_t>(2) ||
          SplineType::kParametricDimensions == static_cast<std::size_t>(3),
      "MFEM export limited to 2D and 3D meshes");

  // Initialize Vertices to identify their vertex ids
  std::vector<typename SplineType::PhysicalPointType_> corner_vertices;
  corner_vertices.reserve(number_of_splines * number_of_element_vertices);

  // Retrieve all vertices
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    const auto global_vertex_id =
        algorithms::HyperCube<SplineType::kParametricDimensions>::
            VertexIdForDegrees(spline_group[i_spline].GetDegrees());
    const auto& spline_ctps = [&]() {
      if constexpr (isRationalBezierSpline_v<SplineType>) {
        return spline_group[i_spline].GetWeightedControlPoints();
      } else {
        return spline_group[i_spline].control_points;
      }
    }();
    for (std::size_t i_point{}; i_point < number_of_element_vertices;
         i_point++) {
      corner_vertices.push_back(spline_ctps[global_vertex_id[i_point]]);
    }
  }

  const auto mfem_information_tuple = utils::algorithms::ExtractMFEMInformation(
      corner_vertices,
      spline_group.MaximumCorner() - spline_group.MinimumCorner());

  // Write Header
  export_file << "MFEM NURBS mesh v1.0\n";

  // Write dimension
  export_file << "\ndimension \n" << SplineType::kParametricDimensions << "\n";

  // Write segment element
  export_file << "\nelements\n";
  export_file << number_of_splines << "\n";

  // Maximum vertex_id for vertex_output
  const auto& vertex_id = std::get<1>(mfem_information_tuple);
  std::size_t maximum_vertex_id{0};

  // Write Elements
  // structure of line:  <element attribute><geometric type><vertex 0>...
  // <vertex m>
  // run through all splines
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    // write the first two entries of the line:
    // geometric type: 3 = square | 5 = cube
    export_file << 1 << " "
                << ((SplineType::kParametricDimensions == 2) ? 3 : 5);
    // run through all vertex inside the current line
    for (std::size_t i_vertex{}; i_vertex < number_of_element_vertices;
         i_vertex++) {
      const std::size_t& c_vertex =
          vertex_id[i_spline * number_of_element_vertices + i_vertex];
      export_file << " " << c_vertex;
      // Update maximum vertex_id, this is needed inside the vertices segment
      maximum_vertex_id =
          (c_vertex > maximum_vertex_id) ? c_vertex : maximum_vertex_id;
    }
    export_file << "\n";
  }

  // Write out boundaries
  // structure of line:  <boundary id><boundary element type><vertex 0>...
  // <vertex m>
  export_file << "\nboundary\n";

  const auto& boundaries = std::get<3>(mfem_information_tuple);
  const std::size_t number_of_boundaries = boundaries.size();
  assert(number_of_boundaries > 0);
  const std::size_t number_of_vertices_per_boundary_element =
      boundaries[0].size();
  // write the number of boundaries
  export_file << number_of_boundaries << "\n";
  for (std::size_t i_boundary{}; i_boundary < number_of_boundaries;
       i_boundary++) {
    // write the first two entries of the line:
    // geometric type: 1 = line | 3 = square
    export_file << 1 << " "
                << ((SplineType::kParametricDimensions == 2) ? 1 : 3);
    for (std::size_t i_vert{}; i_vert < number_of_vertices_per_boundary_element;
         i_vert++) {
      export_file << " " << boundaries[i_boundary][i_vert];
    }
    export_file << "\n";
  }

  // Export Edges (Corner line segments)
  // structure of line:  <knot vector id><vertex start><vertex end>
  export_file << "\n\nedges\n";
  // get the number of edges from the map
  const auto& edge_information = std::get<2>(mfem_information_tuple);
  const std::size_t n_edges = edge_information.size();
  assert(n_edges > 0);
  // write the number of edges
  export_file << n_edges << "\n";
  for (std::size_t i_edge{}; i_edge < n_edges; i_edge++) {
    // write edges
    export_file << edge_information[i_edge][0] << " "
                << edge_information[i_edge][1] << " "
                << edge_information[i_edge][2] << "\n";
  }
  export_file << "\n";

  // Export Number of Vertices
  export_file << "\n\nvertices\n";
  // write the number of vertices
  export_file << maximum_vertex_id + 1 << "\n";

  // Export all Splines as patches and their control points
  export_file << "\n\npatches\n";
  // run through all splines/patches and write the knot vector- and control
  // point segments
  for (std::size_t i_spline{}; i_spline < number_of_splines; i_spline++) {
    export_file << "\nknotvectors"
                << "\n";
    // write the number of knot vectors
    export_file << SplineType::kParametricDimensions << "\n";

    // structure of line <degree><number of control points><xi 0>...<xi m>
    // copy code from format2XMLfile  ;)
    for (std::size_t i_par_dim{}; i_par_dim < SplineType::kParametricDimensions;
         i_par_dim++) {
      // Write All Control Points
      export_file << spline_group[i_spline].GetDegrees()[i_par_dim]
                  << " "
                  // Number of control points in parametric dimension
                  << spline_group[i_spline].GetDegrees()[i_par_dim] + 1 << " ";
      for (std::size_t i_point{};
           i_point < 2 * spline_group[i_spline].GetDegrees()[i_par_dim] + 2;
           i_point++) {
        export_file
            << (i_point < (spline_group[i_spline].GetDegrees()[i_par_dim] + 1)
                    ? "0.0 "
                    : "1.0 ");
      }
      export_file << "\n";
    }

    // write physical dimension
    export_file << "\ndimension"
                << "\n";
    export_file << SplineType::kParametricDimensions << "\n";

    // write the control point segment
    // structure of line <x><y>(<z>)<weight>
    export_file << "\n\ncontrolpoints_cartesian\n";
    for (std::size_t i_point{};
         i_point < spline_group[i_spline].GetNumberOfControlPoints();
         i_point++) {
      for (std::size_t i_dim{};
           i_dim < SplineType::PhysicalPointType_ ::kSpatialDimension;
           i_dim++) {
        const auto local_ctp = [&]() {
          if constexpr (isRationalBezierSpline_v<SplineType>) {
            return spline_group[i_spline].GetWeightedControlPoints()[i_point] /
                   spline_group[i_spline].GetWeights()[i_point];
          } else {
            return spline_group[i_spline].control_points[i_point];
          }
        }();
        export_file << local_ctp[i_dim] << " ";
      }

      // add weights (polynomial splines are always unweighted)
      export_file << [&]() {
        if constexpr (isRationalBezierSpline_v<SplineType>) {
          return spline_group[i_spline].GetWeights()[i_point];
        } else {
          return 1.;
        }
      }() << "\n";
    }

    export_file << "\n";
  }

  // Write Footer
  export_file.close();
  Logger::UserInfo("Finised XML Export");
}

#endif  // SRC_UTILS_EXPORT_INC
