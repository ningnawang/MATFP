// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/Common.h"

#include <geogram/basic/permutation.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_reorder.h>

#include <algorithm>
#include <fstream>
#include <stack>

#include "matfp/Args.h"
#include "matfp/Logger.h"
#include "matfp/MeshIO.h"

namespace matfp {
/*
 * Get File Name from a Path with or without extension
 */

std::string getFileExtension(std::string filePath) {
  return filePath.substr(filePath.find_last_of(".") + 1);
}

std::string getFileName(std::string filePath, bool withExtension,
                        char seperator) {
  std::string filename_ext =
      filePath.substr(filePath.find_last_of(seperator) + 1);
  if (withExtension) return filename_ext;
  size_t lastindex = filename_ext.find_last_of(".");
  return filename_ext.substr(0, lastindex);
}

void sample_triangle(const std::array<Vector3, 3>& vs,
                     std::vector<GEO::vec3>& ps, Scalar sampling_dist) {
  Scalar sqrt3_2 = std::sqrt(3) / 2;

  std::array<Scalar, 3> ls;
  for (int i = 0; i < 3; i++) {
    ls[i] = (vs[i] - vs[mod3(i + 1)]).squaredNorm();
  }
  auto min_max = std::minmax_element(ls.begin(), ls.end());
  int min_i = min_max.first - ls.begin();
  int max_i = min_max.second - ls.begin();
  Scalar N = sqrt(ls[max_i]) / sampling_dist;
  if (N <= 1) {
    for (int i = 0; i < 3; i++)
      ps.push_back(GEO::vec3(vs[i][0], vs[i][1], vs[i][2]));
    return;
  }
  if (N == int(N)) N -= 1;

  GEO::vec3 v0(vs[max_i][0], vs[max_i][1], vs[max_i][2]);
  GEO::vec3 v1(vs[mod3(max_i + 1)][0], vs[mod3(max_i + 1)][1],
               vs[mod3(max_i + 1)][2]);
  GEO::vec3 v2(vs[mod3(max_i + 2)][0], vs[mod3(max_i + 2)][1],
               vs[mod3(max_i + 2)][2]);

  GEO::vec3 n_v0v1 = GEO::normalize(v1 - v0);
  for (int n = 0; n <= N; n++) {
    ps.push_back(v0 + n_v0v1 * sampling_dist * n);
  }
  ps.push_back(v1);

  Scalar h = GEO::distance(
      GEO::dot((v2 - v0), (v1 - v0)) * (v1 - v0) / ls[max_i] + v0, v2);
  int M = h / (sqrt3_2 * sampling_dist);
  if (M < 1) {
    ps.push_back(v2);
    return;
  }

  GEO::vec3 n_v0v2 = GEO::normalize(v2 - v0);
  GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
  Scalar tan_v0, tan_v1, sin_v0, sin_v1;
  sin_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) /
           (GEO::distance(v0, v2) * GEO::distance(v0, v1));
  tan_v0 = GEO::length(GEO::cross((v2 - v0), (v1 - v0))) /
           GEO::dot((v2 - v0), (v1 - v0));
  tan_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) /
           GEO::dot((v2 - v1), (v0 - v1));
  sin_v1 = GEO::length(GEO::cross((v2 - v1), (v0 - v1))) /
           (GEO::distance(v1, v2) * GEO::distance(v0, v1));

  for (int m = 1; m <= M; m++) {
    int n = sqrt3_2 / tan_v0 * m + 0.5;
    int n1 = sqrt3_2 / tan_v0 * m;
    if (m % 2 == 0 && n == n1) {
      n += 1;
    }
    GEO::vec3 v0_m = v0 + m * sqrt3_2 * sampling_dist / sin_v0 * n_v0v2;
    GEO::vec3 v1_m = v1 + m * sqrt3_2 * sampling_dist / sin_v1 * n_v1v2;
    if (GEO::distance(v0_m, v1_m) <= sampling_dist) break;

    Scalar delta_d =
        ((n + (m % 2) / 2.0) - m * sqrt3_2 / tan_v0) * sampling_dist;
    GEO::vec3 v = v0_m + delta_d * n_v0v1;
    int N1 = GEO::distance(v, v1_m) / sampling_dist;
    for (int i = 0; i <= N1; i++) {
      ps.push_back(v + i * n_v0v1 * sampling_dist);
    }
  }
  ps.push_back(v2);

  // sample edges
  N = sqrt(ls[mod3(max_i + 1)]) / sampling_dist;
  if (N > 1) {
    if (N == int(N)) N -= 1;
    GEO::vec3 n_v1v2 = GEO::normalize(v2 - v1);
    for (int n = 1; n <= N; n++) {
      ps.push_back(v1 + n_v1v2 * sampling_dist * n);
    }
  }

  N = sqrt(ls[mod3(max_i + 2)]) / sampling_dist;
  if (N > 1) {
    if (N == int(N)) N -= 1;
    GEO::vec3 n_v2v0 = GEO::normalize(v0 - v2);
    for (int n = 1; n <= N; n++) {
      ps.push_back(v2 + n_v2v0 * sampling_dist * n);
    }
  }
}

int get_t(const Vector3& p0, const Vector3& p1, const Vector3& p2) {
  static const std::array<Vector3, 3> ns = {
      {Vector3(1, 0, 0), Vector3(0, 1, 0), Vector3(0, 0, 1)}};

  Vector3 n = (p1 - p2).cross(p0 - p2);
  Scalar max = 0;
  int t = 0;
  for (int i = 0; i < 3; i++) {
    //        Scalar cos_a = abs(n.dot(ns[i]));
    Scalar cos_a = abs(n[i]);
    if (cos_a > max) {
      max = cos_a;
      t = i;
    }
  }
  return t;
}
Vector2 to_2d(const Vector3& p, int t) {
  return Vector2(p[(t + 1) % 3], p[(t + 2) % 3]);
}
Vector2 to_2d(const Vector3& p, const Vector3& n, const Vector3& pp, int t) {
  Scalar dist = n.dot(p - pp);
  Vector3 proj_p = p - dist * n;
  return Vector2(proj_p[(t + 1) % 3], proj_p[(t + 2) % 3]);
}

Vector3 get_normal(const Vector3& a, const Vector3& b, const Vector3& c) {
  // facet abc is clockwise/counter-clockwise
  // (a-c) x (b-c) is clockwise/counter-clockwise
  return ((a - c).cross(b - c)).normalized();
}

// unit vector from a to b
Vector3 get_direction(const Vector3& a, const Vector3& b) {
  return (b - a).normalized();
}

bool is_vector_same(const Vector3& a, const Vector3& b, const double eps) {
  return (a - b).norm() <= eps;
}

bool is_vector_same_direction(const Vector3& a, const Vector3& b) {
  // angle betwen a and b is in [0,90] or [270,360]
  // => cos(angle) >= 0
  return a.dot(b) / (a.norm() * b.norm()) >= SCALAR_ZERO;
}

bool is_vector_same_direction(const Vector3& a, const Vector3& b,
                              const double degree) {
  // angle betwen a and b is in [0, degree]
  // const double halfC = M_PI / 180;
  return a.dot(b) / (a.norm() * b.norm()) >= std::cos(degree * HALF_PI);
}

bool is_vector_oppo_direction(const Vector3& a, const Vector3& b,
                              const double degree) {
  // angle betwen a and b is in [180-degree, 180]
  double cos_value = a.dot(b) / (a.norm() * b.norm());
  return cos_value <= std::cos((180. - degree) * HALF_PI) && cos_value >= -1;
}

// check if x is within range of (a, b) vectors
bool is_vector_within_range_of_two_vectors(const Vector3& a, const Vector3& b,
                                           const Vector3& x) {
  // https://stackoverflow.com/a/43384516
  Vector3 ab = a.cross(b);  // right system, angle from a to b
  Vector3 ax = a.cross(x);
  Vector3 xb = x.cross(b);
  if (is_vector_same_direction(ab, ax, 90) &&
      is_vector_same_direction(ab, xb, 90))
    return true;
  return false;
}

// Two unit vectors can be considered as two points on unit sphere.
// On a unit sphere, the spherical distance between two points
// is just the radian between them.
// here the spherical distance of two unit vectors d(p,q) = arcos(p.dot(q))
// the domain of the arccos function is from [−1, +1]
// and the range is from [0, pi] radians (or from 0° to 180°).
double spherical_dist_two_unit_vectors(const Vector3& v1, const Vector3& v2) {
  return std::acos(v1.dot(v2));  // always positive
}

double angle_between_two_vectors_in_degrees(const GEO::vec3& a,
                                            const GEO::vec3& b) {
  // a and b should be normalized, so a.dot(b) in range [-1, 1]
  // However, with floating-point math, this is not necessarily true
  // so we need to clamp into [-1, 1]
  GEO::vec3 a_normalized = GEO::normalize(a);
  GEO::vec3 b_normalized = GEO::normalize(b);
  double ab_dot = GEO::dot(a_normalized, b_normalized);
  if (ab_dot <= -1.0) {
    return 180.;
  } else if (ab_dot >= 1.0) {
    return 0.;
  }
  // ab_dot in (-1, 1)
  double diff_angle = std::acos(ab_dot);
  // logger().debug("calculate angle ({},{},{}), ({},{},{}), diff_angle {}",
  //     a[0], a[1], a[2], b[0], b[1], b[2], diff_angle
  // );
  diff_angle *= (180. / PI);
  return diff_angle;
}

double angle_between_two_vectors_in_degrees(const Vector3& a,
                                            const Vector3& b) {
  // a and b should be normalized, so a.dot(b) in range [-1, 1]
  // However, with floating-point math, this is not necessarily true
  // so we need to clamp into [-1, 1]
  Vector3 a_normalized = a.normalized();
  Vector3 b_normalized = b.normalized();
  double ab_dot = a_normalized.dot(b_normalized);
  if (ab_dot <= -1.0) {
    return 180.;
  } else if (ab_dot >= 1.0) {
    return 0.;
  }
  // ab_dot in (-1, 1)
  double diff_angle = std::acos(ab_dot);
  // logger().debug("calculate angle ({},{},{}), ({},{},{}), diff_angle {}",
  //     a[0], a[1], a[2], b[0], b[1], b[2], diff_angle
  // );
  diff_angle *= (180. / PI);
  return diff_angle;
}

double get_distance_between_two_vectors(const Vector3& a, const Vector3& b) {
  return (a - b).norm();
}

double get_distance_to_a_line(const Vector3& p, const Vector3& line_p,
                              const Vector3& line_dir) {
  double dist_2_proj = std::abs((p - line_p).dot(line_dir));
  double dist_2_p = get_distance_between_two_vectors(p, line_p);
  double dist_2_line =
      std::sqrt(std::pow(dist_2_p, 2) - std::pow(dist_2_proj, 2));
  return dist_2_line;
}

// Plane is given by a point p and normal pn,
// return the distance from point q to plane
double get_dist_to_a_plane(const Vector3& p, const Vector3& pn,
                           const Vector3& q) {
  return std::abs((p - q).dot(pn));
}

double sphere_dist_to_plane(const Vector3& center, const double& radius,
                            const Vector3& pl_point, const Vector3& pl_normal) {
  double dist = std::abs(radius - (pl_point - center).dot(pl_normal));
  return dist;
}

double sphere_dist_to_point(const Vector3& center, const double& radius,
                            const Vector3& pl_point, const Vector3& pl_normal) {
  Vector3 proj = center + pl_normal * radius;
  double dist = (proj - pl_point).norm();
  dist = std::abs(dist);
  return dist;
}

Eigen::Matrix4d get_transformation_matrix(const double rangle,
                                          const Vector3& raxis,
                                          const Vector3& translate_vec) {
  // Define our own transformation matrix
  // First rotate, then translate
  double c = std::cos(rangle);
  double s = std::sin(rangle);
  Eigen::Matrix4d m_myself;
  m_myself.row(0) =
      Vector4(raxis[0] * raxis[0] * (1. - c) + c,
              raxis[0] * raxis[1] * (1. - c) - raxis[2] * s,
              raxis[0] * raxis[2] * (1 - c) + raxis[1] * s, translate_vec[0]);
  m_myself.row(1) =
      Vector4(raxis[1] * raxis[0] * (1. - c) + raxis[2] * s,
              raxis[1] * raxis[1] * (1. - c) + c,
              raxis[1] * raxis[2] * (1 - c) - raxis[0] * s, translate_vec[1]);
  m_myself.row(2) =
      Vector4(raxis[2] * raxis[0] * (1. - c) - raxis[1] * s,
              raxis[2] * raxis[1] * (1. - c) + raxis[0] * s,
              raxis[2] * raxis[2] * (1 - c) + c, translate_vec[2]);
  m_myself.row(3) = Vector4(0, 0, 0, 1);
  return m_myself;
}

bool is_eigen_transform_in_sane(
    const Eigen::Transform<double, 3, Eigen::Affine>& trans,
    const double rangle, const Vector3& raxis, const Vector3& translate_vec) {
  // For sanity check
  Eigen::Matrix4d m_eigen = trans.matrix();
  Eigen::Matrix4d m_myself =
      get_transformation_matrix(rangle, raxis, translate_vec);

  logger().debug("--------- Matrix Comparison: ");
  for (int i = 0; i < 3; i++) {
    logger().debug("m_myself.row({}):  ({},{},{},{})", i, m_myself.row(i)[0],
                   m_myself.row(i)[1], m_myself.row(i)[2], m_myself.row(i)[3]);
    logger().debug("m_eigen.row({}):   ({},{},{},{})", i, m_eigen.row(i)[0],
                   m_eigen.row(i)[1], m_eigen.row(i)[2], m_eigen.row(i)[3]);
    if ((m_myself.row(i) - m_myself.row(i)).norm() > SCALAR_ZERO) {
      log_and_throw("ERROR: Matrix is not identical!!!");
      // return false;
    }
  }
  return true;
}

void get_random_vector_between_two_vectors(const Vector3& nA, const Vector3& nB,
                                           Vector3& nX) {
  double rand_k = ((double)std::rand() / (RAND_MAX));
  nX = rand_k * nA + (1. - rand_k) * nB;
  nX.normalize();
  // logger().debug("rand_k: {}, nA ({},{},{}), nB: ({},{},{}), nX: ({},{},{})",
  // 	rand_k, nA[0], nA[1], nA[2], nB[0], nB[1], nB[2], nX[0], nX[1], nX[2]);
}

// k >= 2 includes 2 boundary normals
void sample_k_vectors_given_two_vectors(const Vector3& nA, const Vector3& nB,
                                        const int k,
                                        std::vector<Vector3>& nXs) {
  nXs.clear();
  if (k < 2) return;
  double step = 1. / std::max(k - 1, 1);
  for (double t = 0; t <= 1.; t += step) {
    Vector3 nX = (1. - t) * nA + t * nB;
    nX.normalize();
    nXs.push_back(nX);
  }
  if (nXs.size() == k - 1) {
    double t = 1;
    Vector3 nX = (1. - t) * nA + t * nB;
    nX.normalize();
    nXs.push_back(nX);
  }
  if (k != nXs.size()) {
    logger().error(
        "[ERROR] did not created k: {} normlas with step {}, created {}", k,
        step, nXs.size());
    log_and_throw("ERROR");
  }
}

// add k new sample between segment
// the total samples returned are k+2
void sample_k_between_segment(const Vector3& seg0, const Vector3& seg1,
                              const int k, std::vector<Vector3>& all_samples) {
  all_samples.clear();
  Vector3 vec01 = seg1 - seg0;
  double step = 1. / std::max(1, k + 1);
  for (double t = 0; t <= 1.; t += step) {
    all_samples.push_back(seg0 + t * vec01);
  }
  if (all_samples.size() == k + 2 - 1) {
    // add last one
    all_samples.push_back(seg1);
  }
  if (k + 2 != all_samples.size()) {
    logger().error(
        "[ERROR] did not created k+2: {} samples with step {}, created {}",
        k + 2, step, all_samples.size());
    log_and_throw("ERROR");
  }
}

Vector3 get_centroid(const std::vector<Vector3>& vs) {
  Vector3 c(0., 0., 0.);
  if (vs.size() == 0) return c;
  for (int i = 0; i < vs.size(); i++) {
    c += vs[i];
  }
  c /= vs.size();
  return c;
};

// Return normalized vector
Vector3 get_mesh_facet_normal(const GEO::Mesh& mesh, const int fidx) {
  GEO::vec3 fn = GEO::Geom::mesh_facet_normal(mesh, fidx);
  Vector3 mid_normal = to_eigen(fn);
  mid_normal.normalize();
  return mid_normal;
}

Vector3 get_triangle_centroid(const Vector3& v1, const Vector3& v2,
                              const Vector3& v3) {
  return (v1 + v2 + v3) / 3.;
};

double get_triangle_area(const Vector3& v1, const Vector3& v2,
                         const Vector3& v3) {
  // Heron's formula
  double l1 = (v1 - v2).norm();
  double l2 = (v1 - v3).norm();
  double l3 = (v2 - v3).norm();
  double p = (l1 + l2 + l3) / 2;
  double s = std::sqrt(p * (p - l1) * (p - l2) * (p - l3));
  return s;
}

double get_tetra_volumn(const Vector3& v1, const Vector3& v2, const Vector3& v3,
                        const Vector3& v4) {
  // https://math.stackexchange.com/a/1919737
  // 1/6 * volume of parallelepiped
  // (https://mathinsight.org/image/volume_parallelepiped) volume of
  // parallelpiiped = |(a×b)⋅c|
  Vector3 a = v2 - v1;
  Vector3 b = v3 - v1;
  Vector3 c = v4 - v1;
  double volume = (a.cross(b)).dot(c);
  return 1. / 6. * std::abs(volume);
}

Vector3 to_eigen(const GEO::vec3& v) { return Vector3(v[0], v[1], v[2]); }

GEO::vec3 to_vec3(const Vector3& v) { return GEO::vec3(v[0], v[1], v[2]); }

void print_mesh_info(const GEO::Mesh& M) {
  logger().debug("------ Mesh Info ------");
  logger().debug("cell facets: {}", M.cell_facets.nb());
  logger().debug("cells: {}", M.cells.nb());
  logger().debug("vertex number: {}", M.vertices.nb());
  logger().debug("facet number: {}", M.facets.nb());
  logger().debug("facet corners: {}", M.facet_corners.nb());
}

// Project a given point a on plane.
// Plane is defined by point c and normal n
// return a_projected
// return dist, distance between a and a_projected
void project_point_on_plane(
    const Vector3& c,  // point c and normal n define a plane
    const Vector3& n, const Vector3& a, Vector3& a_projected, double& dist) {
  double t = (c - a).dot(n) / n.dot(n);
  a_projected = a + t * n;
  dist = (a - a_projected).norm();
}

void project_point_onto_line(const Vector3& p, const Vector3& v0,
                             const Vector3& v1, Vector3& p_proj, double& dist) {
  Vector3 v0p(p - v0), v0v1(v1 - v0);
  p_proj = v0 + v0p.dot(v0v1) / v0v1.dot(v0v1) * v0v1;
  dist = (p - p_proj).norm();
}

void color_addon_helper(double& r, double& g, double& b) {
  if (r == 0. && g == 0 && b == 0) {
    // r = 0.01; g = 0.4; b = 0.6;
    r = ((double)rand() / (RAND_MAX));
    g = ((double)rand() / (RAND_MAX));
    b = ((double)rand() / (RAND_MAX));
    return;
  }
  r += 0.1, g += 0.02, b += 0.3;
  if (r >= 1.) r = 0.01;
  if (g >= 1.) g = 0.02;
  if (b >= 1.) b = 0.03;
}

#include <ctime>
#include <string>
std::string get_timestamp() {
  auto now = std::time(nullptr);
  std::ostringstream os;
  // os << std::put_time(std::gmtime(&now),"%F  %T");
  os << std::put_time(std::localtime(&now), "%F_%T");
  return os.str();
}

// rewrite GEO::Geom::tetra_circum_center()
bool tetra_circum_center(const GEO::vec3& p, const GEO::vec3& q,
                         const GEO::vec3& r, const GEO::vec3& s,
                         GEO::vec3& cent_out, bool is_debug) {
  GEO::vec3 qp = q - p;
  double qp2 = GEO::length2(qp);
  GEO::vec3 rp = r - p;
  double rp2 = GEO::length2(rp);
  GEO::vec3 sp = s - p;
  double sp2 = GEO::length2(sp);

  double num_x = GEO::det3x3(qp.y, qp.z, qp2, rp.y, rp.z, rp2, sp.y, sp.z, sp2);
  double num_y = GEO::det3x3(qp.x, qp.z, qp2, rp.x, rp.z, rp2, sp.x, sp.z, sp2);
  double num_z = GEO::det3x3(qp.x, qp.y, qp2, rp.x, rp.y, rp2, sp.x, sp.y, sp2);
  double den =
      GEO::det3x3(qp.x, qp.y, qp.z, rp.x, rp.y, rp.z, sp.x, sp.y, sp.z);

  if (is_debug) {
    logger().debug("tet has den: {}", std::fabs(den));
  }

  // if (std::fabs(den) <= 1e-30) {
  if (std::fabs(den) <= 1e-12) {
    return false;
  }
  den *= 2.0;

  cent_out = GEO::vec3(p.x + num_x / den, p.y - num_y / den, p.z + num_z / den);
  return true;
}

/****************************************************************************/

void set_intersection(std::vector<std::array<double, 2>>& intervals,
                      std::array<double, 2>& range) {
  if (intervals.empty()) return;
  // First interval
  double l = intervals[0][0];
  double r = intervals[0][1];

  // Check rest of the intervals and find the intersection
  for (int i = 1; i < intervals.size(); i++) {
    // If no intersection exists
    if (intervals[i][0] > r || intervals[i][1] < l) {
      logger().error("No intersection found");
      return;
    } else {  // Else update the intersection
      l = std::max(l, intervals[i][0]);
      r = std::min(r, intervals[i][1]);
    }
  }
  range = {{l, r}};
  // logger().debug("intersection range [{},{}]", l, r);
}

bool set_intersection(std::array<double, 2>& interval1,
                      std::array<double, 2>& interval2,
                      std::array<double, 2>& range) {
  std::vector<std::array<double, 2>> intervals;
  intervals.push_back(interval1);
  intervals.push_back(interval2);

  if (intervals.empty()) return false;

  // First interval
  double l = intervals[0][0];
  double r = intervals[0][1];

  // Check rest of the intervals and find the intersection
  for (int i = 1; i < intervals.size(); i++) {
    // If no intersection exists
    // we don't consider equal cases
    if (intervals[i][0] >= r || intervals[i][1] <= l) {
      return false;
    } else {  // Else update the intersection
      l = std::max(l, intervals[i][0]);
      r = std::min(r, intervals[i][1]);
    }
  }
  range = {{l, r}};
  // logger().debug("intersection range [{},{}]", l, r);
  return true;
}

void set_union(std::vector<std::array<double, 2>>& intervals,
               std::vector<std::array<double, 2>>& ranges) {
  if (intervals.empty()) return;
  // sort the intervals in increasing order of start time
  std::sort(intervals.begin(), intervals.end());

  // Create an empty stack of intervals
  std::stack<std::array<double, 2>> s;
  // push the first interval to stack
  s.push(intervals[0]);
  // Start from the next interval and merge if necessary
  for (int i = 1; i < intervals.size(); i++) {
    // get interval from stack top
    std::array<double, 2> top = s.top();
    // if current interval is not overlapping with stack top,
    // push it to the stack
    if (top[1] < intervals[i][0]) s.push(intervals[i]);
    // Otherwise update the ending time of top if ending of current
    // interval is more
    else if (top[1] < intervals[i][1]) {
      top[1] = intervals[i][1];
      s.pop();
      s.push(top);
    }
  }
  // no need to clear ranges, it might contains value
  while (!s.empty()) {
    ranges.push_back(s.top());
    s.pop();
  }
  // logger().debug("union range {}", ranges);
  return;
}

double get_median(std::vector<double> vec) {  // copy here, yes, not reference
  size_t size = vec.size();
  if (size == 0) {
    log_and_throw("ERROR: cannot find median for vector, size = 0");
  }
  std::sort(vec.begin(), vec.end());
  if (size % 2 == 0) {
    return (vec[size / 2 - 1] + vec[size / 2]) / 2;
  } else {
    return vec[size / 2];
  }
}

double get_average(const std::vector<double>& vec) {
  size_t size = vec.size();
  if (size == 0) {
    log_and_throw("ERROR: cannot find average for vector, size = 0");
  }
  double sum = vec[0];
  for (int i = 1; i < vec.size(); i++) {
    sum += vec[i];
  }
  sum /= vec.size();
  return sum;
}

}  // namespace matfp
