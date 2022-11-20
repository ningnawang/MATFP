#pragma once

#include <geogram/basic/geometry.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/mesh/mesh_tetrahedralize.h>

#include <unordered_set>
#include <vector>

#include "matfp/AABBWrapper.h"
#include "matfp/Types/CommonTypes.h"

#define TIMING_BREAKDOWN true

/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */
namespace GEOCommon {

using namespace GEO;

const index_t pid_size = 6;
const int pid[pid_size][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2},
                              {1, 2, 0, 3}, {1, 3, 2, 0}, {2, 3, 0, 1}};

/**
 * \brief Computes the circumcenter and squared radius of a tetrahedron.
 * \param[in] p first vertex of the tetrahedron
 * \param[in] q second vertex of the tetrahedron
 * \param[in] r third vertex of the tetrahedron
 * \param[in] s fourth vertex of the tetrahedron
 * \param[out] t_circumcenter computed circumcenter of the tetrahedron
 * \param[out] squared_radius computed squared radius of the tetrahedron
 * \param[out] dihedral_angles if non-nullptr,
 *   computed dihedral angles of the tetrahedron
 * \return true if successful, false otherwise
 *  (for instance, if the tetrahedron is flat).
 */
inline bool tetra_circumcenter_squaredradius(
    const vec3& p, const vec3& q, const vec3& r, const vec3& s,
    vec3& t_circumcenter, double& squared_radius,
    double* dihedral_angles = nullptr) {
  vec3 qp = q - p;
  vec3 rp = r - p;
  vec3 sp = s - p;

  double tet_vol_x_12 = 2.0 * dot(qp, cross(rp, sp));

  // it is not a safe check, one should avoid
  // to input a degenerated tetra.
  if (tet_vol_x_12 == 0.0) {
    return false;
  }

  vec3 d_t_circumcenter = length2(qp) * cross(rp, sp) +
                          length2(rp) * cross(sp, qp) +
                          length2(sp) * cross(qp, rp);

  squared_radius = length2(d_t_circumcenter) / (tet_vol_x_12 * tet_vol_x_12);
  d_t_circumcenter = (1.0 / tet_vol_x_12) * d_t_circumcenter;
  t_circumcenter = p + d_t_circumcenter;

  if (dihedral_angles) {
    // compute Dihedral angle directly
    // the following computation is not optimized.
    vec3 point[4];
    point[0] = p;
    point[1] = q;
    point[2] = r;
    point[3] = s;
    for (index_t i = 0; i < pid_size; i++) {
      vec3 b1 = point[pid[i][0]] - point[pid[i][2]];
      vec3 b2 = point[pid[i][1]] - point[pid[i][0]];
      vec3 b3 = point[pid[i][3]] - point[pid[i][1]];
      vec3 b2b3 = cross(b2, b3);
      dihedral_angles[i] =
          ::fabs(::atan2(length(b2) * dot(b1, b2b3), dot(cross(b1, b2), b2b3)));
    }
  }
  return true;
}

/**
 * \brief Gets a Delaunay vertex by global index.
 * \param[in] delaunay the Delaunay triangulation
 * \param[in] v the index of the vertex
 * \return a const reference to the vertex, as a vec3
 */
inline const vec3& delaunay_vertex(Delaunay* delaunay, index_t v) {
  return *(const vec3*)delaunay->vertex_ptr(v);
}

/**
 * \brief Gets a Delaunay vertex by tetrahedron index
 *  and local vertex index.
 * \param[in] delaunay the Delaunay triangulation
 * \param[in] c the index of the tetrahedron
 * \param[in] lv the local index of the vertex (0,1,2 or 3)
 *  in tetrahedron \p c.
 * \return a const reference to the vertex, as a vec3.
 */
inline const vec3& delaunay_tet_vertex(Delaunay* delaunay, index_t c,
                                       index_t lv) {
  return delaunay_vertex(delaunay, index_t(delaunay->cell_vertex(c, lv)));
}

const index_t tet_facet_vertex[4][3] = {
    {1, 2, 3}, {0, 3, 2}, {3, 0, 1}, {2, 1, 0}};

/**
 * \brief Computes the normal to a tetrahedron facet.
 * \param[in] delaunay the Delaunay triangulation
 * \param[in] t the index of the tetrahedron
 * \param[in] f the local index (0,1,2 or 3) of
 *  the facet in the tetrahedron \p t
 * \return the normal to the facet \p f or tetrahedron \p t
 */
inline vec3 delaunay_facet_normal(Delaunay* delaunay, index_t t, index_t f) {
  geo_debug_assert(f < 4);
  index_t v1 = index_t(delaunay->cell_vertex(t, tet_facet_vertex[f][0]));
  index_t v2 = index_t(delaunay->cell_vertex(t, tet_facet_vertex[f][1]));
  index_t v3 = index_t(delaunay->cell_vertex(t, tet_facet_vertex[f][2]));
  const vec3& p1 = delaunay_vertex(delaunay, v1);
  const vec3& p2 = delaunay_vertex(delaunay, v2);
  const vec3& p3 = delaunay_vertex(delaunay, v3);
  return cross(p2 - p1, p3 - p1);
}

}  // namespace GEOCommon

namespace matfp {

/////////////// I/O
/*
 * Get File Name from a Path with or without extension
 */
std::string getFileExtension(std::string filePath);
std::string getFileName(std::string filePath, bool withExtension = true,
                        char seperator = '/');

/////////////// matrix
inline void remove_row(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
  unsigned int numRows = matrix.rows() - 1;
  unsigned int numCols = matrix.cols();
  if (rowToRemove < numRows)
    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) =
        matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);
  matrix.conservativeResize(numRows, numCols);
}

inline void remove_column(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols() - 1;
  if (colToRemove < numCols)
    matrix.block(0, colToRemove, numRows, numCols - colToRemove) =
        matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);
  matrix.conservativeResize(numRows, numCols);
}

///////////////
void sample_triangle(const std::array<Vector3, 3>& vs,
                     std::vector<GEO::vec3>& ps, Scalar sampling_dist);
bool sample_triangle_and_check_is_out(const std::array<Vector3, 3>& vs,
                                      Scalar sampling_dist, Scalar eps_2,
                                      const AABBWrapper& tree,
                                      GEO::index_t& prev_facet);

////////// Surface
Vector3 get_centroid(const std::vector<Vector3>& vs);
Vector3 get_mesh_facet_normal(const GEO::Mesh& mesh, const int fidx);
Vector3 get_triangle_centroid(const Vector3& v1, const Vector3& v2,
                              const Vector3& v3);
double get_triangle_area(const Vector3& v1, const Vector3& v2,
                         const Vector3& v3);
double get_tetra_volumn(const Vector3& v1, const Vector3& v2, const Vector3& v3,
                        const Vector3& v4);

Vector3 get_normal(const Vector3& a, const Vector3& b, const Vector3& c);
Vector3 get_direction(const Vector3& a, const Vector3& b);
bool is_vector_same(const Vector3& a, const Vector3& b,
                    const double eps = SCALAR_ZERO);
bool is_vector_same_direction(const Vector3& a, const Vector3& b);
bool is_vector_same_direction(const Vector3& a, const Vector3& b,
                              const double degree);
bool is_vector_oppo_direction(const Vector3& a, const Vector3& b,
                              const double degree);
bool is_vector_within_range_of_two_vectors(const Vector3& a, const Vector3& b,
                                           const Vector3& x);
double spherical_dist_two_unit_vectors(const Vector3& v1, const Vector3& v2);
double angle_between_two_vectors_in_degrees(const GEO::vec3& a,
                                            const GEO::vec3& b);
double angle_between_two_vectors_in_degrees(const Vector3& a, const Vector3& b);
double get_distance_between_two_vectors(const Vector3& a, const Vector3& b);
double get_distance_to_a_line(const Vector3& p, const Vector3& line_p,
                              const Vector3& line_dir);
double get_dist_to_a_plane(const Vector3& p, const Vector3& pn,
                           const Vector3& q);
double sphere_dist_to_plane(const Vector3& center, const double& radius,
                            const Vector3& pl_point, const Vector3& pl_normal);
double sphere_dist_to_point(const Vector3& center, const double& radius,
                            const Vector3& pl_point, const Vector3& pl_normal);
Eigen::Matrix4d get_transformation_matrix(const double rangle,
                                          const Vector3& raxis,
                                          const Vector3& translate_vec);
bool is_eigen_transform_in_sane(
    const Eigen::Transform<double, 3, Eigen::Affine>& trans,
    const double rangle, const Vector3& raxis, const Vector3& translate_vec);

void get_random_vector_between_two_vectors(const Vector3& nA, const Vector3& nB,
                                           Vector3& nX);

// totol k
void sample_k_vectors_given_two_vectors(const Vector3& nA, const Vector3& nB,
                                        const int k, std::vector<Vector3>& nXs);
// total k+2
void sample_k_between_segment(const Vector3& seg0, const Vector3& seg1,
                              const int k, std::vector<Vector3>& all_samples);
int get_t(const Vector3& p0, const Vector3& p1, const Vector3& p2);
Vector2 to_2d(const Vector3& p, int t);
Vector2 to_2d(const Vector3& p, const Vector3& n, const Vector3& pp, int t);

Vector3 to_eigen(const GEO::vec3& v);
GEO::vec3 to_vec3(const Vector3& v);

void print_mesh_info(const GEO::Mesh& M);

void color_addon_helper(double& r, double& g, double& b);
std::string get_timestamp();

bool tetra_circum_center(const GEO::vec3& p, const GEO::vec3& q,
                         const GEO::vec3& r, const GEO::vec3& s,
                         GEO::vec3& cent_out, bool is_debug = false);

//////////
// Project a point a on a given plane
// 1. n parallel to vec(ca)
// 2. vec(ca) perpendicular to n
void project_point_on_plane(
    const Vector3& c,  // point c and normal n define a plane
    const Vector3& n, const Vector3& a, Vector3& a_projected, double& dist);
//////////
// Project a point a on a given line
void project_point_onto_line(const Vector3& p, const Vector3& v0,
                             const Vector3& v1, Vector3& p_proj, double& dist);

inline int mod3(int j) { return j % 3; }
inline int mod4(int j) { return j % 4; }
inline int modk(int j, int k) { return j % k; }
//////////
// Not been used
// An alternative for atan2(): https://stackoverflow.com/a/16561333
// atan2 has a discontinuity, there's a point where atan2 will jump between -pi
// and +pi. Input:  dx, dy: coordinates of a (difference) vector. Output: a
// number from the range [-2 .. 2] which is monotonic
//         in the angle this vector makes against the x axis.
//         and with the same discontinuity as atan2
inline double pseudoangle(double dx, double dy) {
  return copysign(1. - dx / (fabs(dx) + fabs(dy)), dy);
}
// Driver function to sort the vector elements by
// first element of pair in descending order
inline bool sort_in_rev(const std::pair<double, int>& a,
                        const std::pair<double, int>& b) {
  return (a.first > b.first);
}

// Combination
// from total_vec to choose k elements
template <typename T>
void get_combination(const std::vector<T>& vec, int offset, int k,
                     std::vector<T>& comb,  // tmp to store one combination
                     std::vector<std::vector<T>>& combinations) {
  if (k >= vec.size()) {
    logger().error("ERROR: k: {}, n: {}", k, vec.size());
    log_and_throw("Error: combination k >= n");
  }

  if (k == 0) {
    combinations.push_back(comb);
    return;
  }
  for (int i = offset; i <= vec.size() - k; i++) {
    comb.push_back(vec[i]);
    get_combination(vec, i + 1, k - 1, comb, combinations);
    comb.pop_back();
  }
}

template <typename T>
void vector_unique(std::vector<T>& v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

double get_median(std::vector<double> vec);
double get_average(const std::vector<double>& vec);

///////////////
// Functions to get the intersection
//
// NOTE:
// The compiler must be able to see the implementation
// in order to generate code for all specializations in your code
// So all template functions must be implemented IN the header

// template <typename T>
// void set_intersection(const std::unordered_set<T>& s1, const
// std::unordered_set<T>& s2, std::vector<T>& v); template <typename T> void
// set_intersection(const std::unordered_set<T>& s1, const
// std::unordered_set<T>& s2, std::unordered_set<T>& v); template <typename T>
// void set_intersection(const std::unordered_set<T>& s1, const
// std::unordered_set<T>& s2, const std::unordered_set<T>& s3, std::vector<T>&
// v); template <typename T> void set_intersection(const std::set<T>& s1, const
// std::set<T>& s2, const std::set<T>& s3, std::vector<T>& v);

// template <typename T>
// void set_intersection(const std::vector<T>& s1, const std::vector<T>& s2,
// std::vector<T>& v); template <typename T> void set_intersection(const
// std::vector<T>& s1, const std::vector<T>& s2, const std::vector<T>& s3,
// std::vector<T>& v); template <typename T> void set_intersection_sorted(const
// std::vector<T>& s1, const std::vector<T>& s2, const std::vector<T>& s3,
// std::vector<T>& v);

template <typename T>
void set_difference(const std::set<T>& s1, const std::set<T>& s2,
                    std::set<T>& v) {
  v.clear();
  for (int x1 : s1) {
    if (s2.find(x1) == s2.end()) {
      v.insert(x1);
    }
  }
  for (int x2 : s2) {
    if (s1.find(x2) == s1.end()) {
      v.insert(x2);
    }
  }
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2, std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2,
                      std::unordered_set<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.insert(x);
    }
  }
}

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2,
                      const std::unordered_set<T>& s3, std::vector<T>& v) {
  if (s2.size() < s1.size() && s2.size() < s1.size()) {
    set_intersection(s2, s1, s3, v);
    return;
  }

  if (s3.size() < s1.size() && s3.size() < s2.size()) {
    set_intersection(s3, s1, s2, v);
    return;
  }

  assert(s1.size() <= s2.size());
  assert(s1.size() <= s3.size());

  v.clear();
  v.reserve(s1.size());
  for (int x : s1) {
    if (s2.count(x) && s3.count(x)) {
      v.push_back(x);
      if (v.size() == 2) break;
    }
  }
}

template <typename T>
void set_intersection(const std::set<T>& s1, const std::set<T>& s2,
                      const std::set<T>& s3, std::vector<T>& v) {
  if (s2.size() < s1.size() && s2.size() < s1.size()) {
    set_intersection(s2, s1, s3, v);
    return;
  }

  if (s3.size() < s1.size() && s3.size() < s2.size()) {
    set_intersection(s3, s1, s2, v);
    return;
  }

  assert(s1.size() <= s2.size());
  assert(s1.size() <= s3.size());

  v.clear();
  v.reserve(s1.size());
  for (int x : s1) {
    if (s2.count(x) && s3.count(x)) {
      v.push_back(x);
      if (v.size() == 2) break;
    }
  }
}

template <typename T>
void set_intersection(const std::vector<T>& s11, const std::vector<T>& s22,
                      std::vector<T>& v) {
  v.clear();
  std::vector<T> s1 = s11;
  std::vector<T> s2 = s22;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(v));
}

template <typename T>
void set_intersection(const std::vector<T>& s11, const std::vector<T>& s22,
                      const std::vector<T>& s33, std::vector<T>& v) {
  v.clear();
  std::vector<T> s1 = s11;
  std::vector<T> s2 = s22;
  std::vector<T> s3 = s33;
  std::sort(s1.begin(), s1.end());
  std::sort(s2.begin(), s2.end());
  std::sort(s3.begin(), s3.end());
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(v));
  auto it = std::set_intersection(v.begin(), v.end(), s3.begin(), s3.end(),
                                  v.begin());
  v.resize(it - v.begin());
}

template <typename T>
void set_intersection_sorted(const std::vector<T>& s1, const std::vector<T>& s2,
                             const std::vector<T>& s3, std::vector<T>& v) {
  std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(v));
  auto it = std::set_intersection(v.begin(), v.end(), s3.begin(), s3.end(),
                                  v.begin());
  v.resize(it - v.begin());
}

void set_intersection(std::vector<std::array<double, 2>>& intervals,
                      std::array<double, 2>& range);
bool set_intersection(std::array<double, 2>& interval1,
                      std::array<double, 2>& interval2,
                      std::array<double, 2>& range);

///////////////
// Functions to get the union
void set_union(std::vector<std::array<double, 2>>& intervals,
               std::vector<std::array<double, 2>>& ranges);

///////////
// Functions to merge two maps
template <typename K, typename V>
std::map<K, std::unordered_set<V>> merge_two_maps(
    std::map<K, std::unordered_set<V>> const& lhs,
    std::map<K, std::unordered_set<V>> const& rhs) {
  typedef typename std::map<K, std::unordered_set<V>>::const_iterator
      input_iterator;
  std::map<K, std::unordered_set<V>> result;
  for (input_iterator it1 = lhs.begin(), it2 = rhs.begin();
       it1 != lhs.end() && it2 != rhs.end();) {
    if (it1->first == it2->first) {
      set_intersection(it1->second, it2->second, result[it1->first]);
      ++it1;
      ++it2;
    } else {
      if (it1->first < it2->first)
        ++it1;
      else
        ++it2;
    }
  }
  return result;
}

template <typename K, typename V>
std::map<K, std::vector<V>> merge_two_maps_give_vectors(
    std::map<K, std::unordered_set<V>> const& lhs,
    std::map<K, std::unordered_set<V>> const& rhs) {
  typedef typename std::map<K, std::unordered_set<V>>::const_iterator
      input_iterator;
  std::map<K, std::vector<V>> result;
  for (input_iterator it1 = lhs.begin(), it2 = rhs.begin();
       it1 != lhs.end() && it2 != rhs.end();) {
    if (it1->first == it2->first) {
      set_intersection(it1->second, it2->second, result[it1->first]);
      ++it1;
      ++it2;
    } else {
      if (it1->first < it2->first)
        ++it1;
      else
        ++it2;
    }
  }
  return result;
}

}  // namespace matfp
