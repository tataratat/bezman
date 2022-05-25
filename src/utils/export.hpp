#ifndef SRC_UTILS_EXPORT_HPP
#define SRC_UTILS_EXPORT_HPP

#include <fstream>
#include <map>
#include <stdexcept>
#include <string>

#include "bezierManipulation/src/bezier_spline.hpp"
#include "bezierManipulation/src/bezier_spline_group.hpp"
#include "bezierManipulation/src/point.hpp"
#include "bezierManipulation/src/utils/algorithms/int_power.hpp"
#include "bezierManipulation/src/utils/algorithms/point_uniquifier.hpp"
#include "bezierManipulation/src/utils/base64.hpp"
#include "bezierManipulation/src/utils/logger.hpp"

namespace beziermanipulation::utils {

/*
 * Static Class for Export Routines
 *
 * Export Routines are extended by whatever format of splines is required.
 *
 * Currently supported:
 *  *.xml -> standard CATS export file (FEAFA, XNS)
 *  *.itd -> IRIT export file extension
 *  *.mesh -> MFEM/GLVis export file
 *  *.json -> Python readable json format (using gustaf Keywords). Also supports
 *            base64 encoding for exact export of the control points
 */
class Export {
 private:
  /*
   * Returns a string with a specific file extension
   *
   * test.txt and test both return test.txt
   */
  static std::string ensureFileExtension(const std::string &filename,
                                         const std::string &obj_extension);

  /*
   * Formats spline to IRIT format
   *
   * Actual implementation of the spline Export to be called for every spline in
   * a group Pipes everything directly into a file
   */
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void format2IRITfile(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
          &spline,
      std::ofstream &export_file);

  /*
   * Formats spline to XML format
   *
   * Actual implementation of the spline Export to be called for every spline in
   * a group Pipes everything directly into a file
   */
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void format2XMLfile(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
          &spline,
      std::ofstream &export_file);

  /**
   * @brief  Formats spline in custom json format
   *
   * Actual implementation of the spline Export to be called for every spline in
   * a group Pipes everything directly into a file
   */
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void format2JSONfile(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
          &spline,
      std::ofstream &export_file, const bool base64encoding);

 public:
  /// Permit creation of class instance
  Export() = delete;

  /// Guess by extension
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void GuessByExtension(
      const BezierSpline<parametric_dimension, PhysicalPointType, ScalarType>
          &spline,
      const std::string &filename);

  /// Guess by extension
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void GuessByExtension(
      const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType> &spline_group,
      const std::string &filename);

  /// Export as IRIT
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsIRIT(const BezierSpline<parametric_dimension, PhysicalPointType,
                                        ScalarType> &spline,
                     const std::string &filename);

  /// Export as IRIT
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsIRIT(
      const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType> &spline_group,
      const std::string &filename);

  /// Export Single Spline as XML
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsXML(const BezierSpline<parametric_dimension, PhysicalPointType,
                                       ScalarType> &spline,
                    const std::string &filename);

  /// Export Group as XML
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsXML(
      const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType> &spline_group,
      const std::string &filename);

  /// Export Single Spline as custom JSON
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsJSON(const BezierSpline<parametric_dimension, PhysicalPointType,
                                        ScalarType> &spline,
                     const std::string &filename,
                     const bool base64encoding = true);

  /// Export Group as  custom JSON
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsJSON(
      const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType> &spline_group,
      const std::string &filename, const bool base64encoding = true);

  /// Export Single Spline MFEM multipatch format
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsMFEM(const BezierSpline<parametric_dimension, PhysicalPointType,
                                        ScalarType> &spline,
                     const std::string &filename);

  /// Export Spline-Group in MFEM multipatch format
  template <std::size_t parametric_dimension, typename PhysicalPointType,
            typename ScalarType>
  static void AsMFEM(
      const BezierSplineGroup<parametric_dimension, PhysicalPointType,
                              ScalarType> &spline_group,
      const std::string &filename);
};

#include "bezierManipulation/src/utils/export.inc"

}  // namespace beziermanipulation::utils
#endif  // SRC_UTILS_EXPORT_HPP
