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

#ifndef SRC_UTILS_EXPORT_HPP
#define SRC_UTILS_EXPORT_HPP

#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>

#include "bezman/src/bezier_group.hpp"
#include "bezman/src/bezier_spline.hpp"
#include "bezman/src/point.hpp"
#include "bezman/src/utils/algorithms/int_power.hpp"
#include "bezman/src/utils/algorithms/point_uniquifier.hpp"
#include "bezman/src/utils/base64.hpp"
#include "bezman/src/utils/logger.hpp"
#include "bezman/src/utils/type_traits/is_bezier_spline.hpp"
#include "bezman/src/utils/type_traits/is_rational_bezier_spline.hpp"

namespace bezman::utils {

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
  template <typename SplineType>
  static void format2IRITfile(const SplineType &spline,
                              std::ofstream &export_file);

  /*
   * Formats spline to XML format
   *
   * Actual implementation of the spline Export to be called for every spline in
   * a group Pipes everything directly into a file
   */
  template <typename SplineType>
  static void format2XMLfile(const SplineType &spline,
                             std::ofstream &export_file);

  /**
   * @brief  Formats spline in custom json format
   *
   * Actual implementation of the spline Export to be called for every spline in
   * a group Pipes everything directly into a file
   */
  template <typename SplineType>
  static void format2JSONfile(const SplineType &spline,
                              std::ofstream &export_file,
                              const bool base64encoding);

 public:
  /// Permit creation of class instance
  Export() = delete;

  /// Guess by extension
  template <typename SplineType>
  static void GuessByExtension(const SplineType &spline,
                               const std::string &filename);

  /// Guess by extension
  template <typename SplineType>
  static void GuessByExtension(const BezierGroup<SplineType> &spline_group,
                               const std::string &filename);

  /// Export as IRIT
  template <typename SplineType>
  static void AsIRIT(const SplineType &spline, const std::string &filename);

  /// Export as IRIT
  template <typename SplineType>
  static void AsIRIT(const BezierGroup<SplineType> &spline_group,
                     const std::string &filename);

  /// Export Single Spline as XML
  template <typename SplineType>
  static void AsXML(const BezierGroup<SplineType> &spline,
                    const std::string &filename);

  /// Export Group as XML
  template <typename SplineType>
  static void AsXML(const SplineType &spline_group,
                    const std::string &filename);

  /// Export Single Spline as custom JSON
  template <typename SplineType>
  static void AsJSON(const SplineType &spline, const std::string &filename,
                     const bool base64encoding = true);

  /// Export Group as  custom JSON
  template <typename SplineType>
  static void AsJSON(const BezierGroup<SplineType> &spline_group,
                     const std::string &filename,
                     const bool base64encoding = true);

  /// Export Single Spline MFEM multipatch format
  template <typename SplineType>
  static void AsMFEM(const SplineType &spline, const std::string &filename);

  /// Export Spline-Group in MFEM multipatch format
  template <typename SplineType>
  static void AsMFEM(const BezierGroup<SplineType> &spline_group,
                     const std::string &filename);
};

#include "bezman/src/utils/export.inc"

}  // namespace bezman::utils
#endif  // SRC_UTILS_EXPORT_HPP
