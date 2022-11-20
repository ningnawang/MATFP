// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "OtherTypes.h"

namespace matfp {
/////////////////////////////////////////////////////////////////////////////////////
// TangentConcaveLine
/////////////////////////////////////////////////////////////////////////////////////
TangentConcaveLine::TangentConcaveLine(const int _id,
                                       const std::array<Vector3, 2> _end_pts,
                                       const std::array<Vector3, 2>& _ns,
                                       const std::array<int, 2>& _fpair,
                                       const std::array<int, 2>& _vpair) {
  id = _id;
  ref_vs_pos = _end_pts;
  direction = (ref_vs_pos[1] - ref_vs_pos[0]).normalized();
  tan_point = (ref_vs_pos[0] + ref_vs_pos[1]) / 2.;
  is_tan_point_updated = false;
  adj_ref_normals = _ns;
  get_random_vector_between_two_vectors(adj_ref_normals[0], adj_ref_normals[1],
                                        normal);
  adj_ref_fs_pair = _fpair;
  ref_vs_seed_ids = _vpair;
  energy = DBL_MAX;
}

bool TangentConcaveLine::operator==(TangentConcaveLine const& b) const {
  Vector3 vec = (b.tan_point - tan_point).normalized();
  if (adj_ref_fs_pair == b.adj_ref_fs_pair) return true;

  // more likely to be different
  // adj_ref_fs_pair might be different, but represent the same line
  if (!is_vector_same_direction(direction, b.direction, esp_degree_5) &&
      !is_vector_oppo_direction(direction, b.direction, esp_degree_5))
    return false;

  if (!is_vector_same_direction(vec, direction, esp_degree_5) &&
      !is_vector_oppo_direction(vec, direction, esp_degree_5))
    return false;

  if (!is_vector_same_direction(vec, b.direction, esp_degree_5) &&
      !is_vector_oppo_direction(vec, b.direction, esp_degree_5))
    return false;

  // vec should parallel to both direction and b.direction
  return true;
}

void TangentConcaveLine::print_info() const {
  logger().debug(
      "---- TangentConcaveLine: energy: {}, energy_over_sq_radius: {}, "
      "ref_vs_seed_ids {}, adj_ref_fs_pair: {}, tan_point ({},{},{}), normal: "
      "({},{},{})",
      energy, energy_over_sq_radius, ref_vs_seed_ids, adj_ref_fs_pair,
      tan_point[0], tan_point[1], tan_point[2], normal[0], normal[1],
      normal[2]);
  logger().debug("adj_ref_normals ({},{},{}) - ({},{},{})",
                 adj_ref_normals[0][0], adj_ref_normals[0][1],
                 adj_ref_normals[0][2], adj_ref_normals[1][0],
                 adj_ref_normals[1][1], adj_ref_normals[1][2]);
}

// if concave line is a curve
// then each line should cover more adjacent faces
bool TangentConcaveLine::is_normal_covered_by_adj_fs(
    const Vector3& n, double esp_degree_given) const {
  if (is_vector_same_direction(adj_ref_normals[0], n, esp_degree_given))
    return true;
  if (is_vector_same_direction(adj_ref_normals[1], n, esp_degree_given))
    return true;
  return false;
}

// given two points of a segment
void TangentConcaveLine::update_tangent_point(const GEO::vec3& v1_p,
                                              const GEO::vec3& v2_p) {
  GEO::vec3 v_mid = 1. / 2. * (v1_p + v2_p);
  Vector3 p2(v_mid[0], v_mid[1], v_mid[2]);
  if (!is_tan_point_updated)
    tan_point = p2;
  else
    tan_point = 1. / 2. * (p2 + tan_point);
}

// static function
double TangentConcaveLine::get_energy_value(const Vector3& theta,
                                            const double& radius,
                                            const double alpha3,
                                            const Vector3& tan_point,
                                            const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  return alpha3 * T.dot(T);
}
// sphere (theta, radiu)
double TangentConcaveLine::update_energy_value(const Vector3& theta,
                                               const double& radius,
                                               const double alpha3) {
  // sanity_check();
  energy = get_energy_value(theta, radius, alpha3, tan_point, normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

/////////////////////////////////////////////////////////////////////////////////////
// TangentPlane
/////////////////////////////////////////////////////////////////////////////////////
TangentPlane::TangentPlane() {
  is_deleted = false;
  is_touched = false;
  extrude_ratio = -1;
  extrude_ratio_signed = 0.;
  fid = -1;
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
};

TangentPlane::TangentPlane(const Vector3& _normal) {
  normal = _normal;
  is_deleted = false;
  is_touched = false;
  extrude_ratio = -1;
  extrude_ratio_signed = 0.;
  fid = -1;
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
}

TangentPlane::TangentPlane(const Vector3& _normal, const Vector3& _point) {
  normal = _normal;
  points.push_back(_point);
  is_deleted = false;
  is_touched = false;
  extrude_ratio = -1;
  extrude_ratio_signed = 0.;
  fid = -1;
  energy = DBL_MAX;
  energy_over_sq_radius = DBL_MAX;
}

bool TangentPlane::operator==(TangentPlane const& b) const {
  return is_vector_same_direction(normal, b.normal, esp_degree_5);
}
/**
 * Compares two 4D points with respect to the lexicographic order.
 * s1, s2 are the coordinates of the two 4D points.
 * true if s1 is strictly before s2 in the lexicographic order.
 * false otherwise.
 */
bool TangentPlane::operator<(const TangentPlane& t2) const {
  if (normal[0] < t2.normal[0]) {
    return true;
  }
  if (normal[0] > t2.normal[0]) {
    return false;
  }
  if (normal[1] < t2.normal[1]) {
    return true;
  }
  if (normal[1] > t2.normal[1]) {
    return false;
  }
  if (normal[2] < t2.normal[2]) {
    return true;
  }
  if (normal[2] > t2.normal[2]) {
    return false;
  }
  return normal[3] < t2.normal[3];
}

bool TangentPlane::clear() {
  points.clear();
  is_deleted = false;
  is_touched = false;
  extrude_ratio = -1;
  fid = -1;
  energy = DBL_MAX;
}

bool TangentPlane::is_same_normal(const Vector3& bnormal,
                                  const double eps_degree) const {
  return is_vector_same_direction(normal, bnormal, eps_degree);
}

bool TangentPlane::is_same_normal(const Vector3& anormal,
                                  const Vector3& bnormal,
                                  const double eps_degree) {
  return is_vector_same_direction(anormal, bnormal, eps_degree);
}

void TangentPlane::union_two_tan_pls_vectors(
    const std::vector<TangentPlane>& tan_pls1,
    const std::vector<TangentPlane>& tan_pls2,
    std::vector<TangentPlane>& tan_pls_merged, const double eps_degree) {
  tan_pls_merged = tan_pls1;
  // check tan_pls2, merge different tangent planes
  for (const auto& tan_pl2 : tan_pls2) {
    bool is_covered = false;
    for (const auto& tan_pl1 : tan_pls_merged) {
      if (tan_pl1.is_same_normal(tan_pl2.normal, eps_degree)) {
        is_covered = true;
        break;
      }
    }  // for tan_pls_merged
    if (!is_covered) tan_pls_merged.push_back(tan_pl2);
  }  // for tan_pls2
}

void TangentPlane::diff_two_tan_pls_vectors(
    const std::vector<TangentPlane>& tan_pls1,
    const std::vector<TangentPlane>& tan_pls2,
    std::vector<TangentPlane>& tan_pls_diff, const double eps_degree) {
  tan_pls_diff.clear();
  // save tan_pl1 if not covered by tan_pls2
  for (const auto& tan_pl1 : tan_pls1) {
    bool is_covered = false;
    for (const auto& tan_pl2 : tan_pls2) {
      if (tan_pl2.is_same_normal(tan_pl1.normal, eps_degree)) {
        is_covered = true;
        break;
      }
    }  // for tan_pls_diff
    if (!is_covered) tan_pls_diff.push_back(tan_pl1);
  }  // for tan_pls1
  // save tan_pl2 if not covered by tan_pls1
  // (NOTE: not tan_pls_diff!)
  for (const auto& tan_pl2 : tan_pls2) {
    bool is_covered = false;
    for (const auto& tan_pl1 : tan_pls1) {
      if (tan_pl1.is_same_normal(tan_pl2.normal, eps_degree)) {
        is_covered = true;
        break;
      }
    }  // for tan_pls1
    if (!is_covered) tan_pls_diff.push_back(tan_pl2);
  }  // for tan_pls2
}

// 1. normals are simialr
// 2. tangent point is close
bool TangentPlane::is_same_tan_pl(const TangentPlane& tan_pl2) const {
  sanity_check();
  tan_pl2.sanity_check();
  if (!is_same_normal(tan_pl2.normal, esp_degree_10)) {
    return false;
  }
  double dist = get_distance_between_two_vectors(points[0], tan_pl2.points[0]);
  if (dist > SCALAR_ZERO_N3) {
    return false;
  }
  return true;
}

void TangentPlane::push_new_point(Vector3 _p) {
  points.push_back(_p);
  is_touched = true;
}

void TangentPlane::print_info() const {
  logger().debug(
      "TangentPlane: is_deleted {}, fid: {} energy: {}, energy_over_sq_radius: "
      "{}, extrude_ratio_signed: {}, points[0] ({},{},{}), normal: "
      "({},{},{}) , points size {}",
      is_deleted, fid, energy, energy_over_sq_radius, extrude_ratio_signed,
      points[0][0], points[0][1], points[0][2], normal[0], normal[1], normal[2],
      points.size());
}

void TangentPlane::update_points_with_centroid() {
  if (points.size() < 2)  // do nothing
    return;
  Vector3 p = get_centroid(points);
  points.clear();
  push_new_point(p);
}

// mainly for updating fids
// but also update points and normal
// make sure they are projected on the surface
void TangentPlane::update_by_sf_mesh(const GEO::Mesh& sf_mesh,
                                     const AABBWrapper& aabb_wrapper) {
  Vector3 near_p;
  update_points_with_centroid();
  int p_fid = aabb_wrapper.project_to_sf_get_nearest_face(points[0]);
  const Vector3 pnormal = get_mesh_facet_normal(sf_mesh, p_fid);
  fid = p_fid;
  normal = pnormal;
}

// static function
double TangentPlane::get_energy_value(const Vector3& theta,
                                      const double& radius, const double alpha1,
                                      const double alpha2,
                                      const Vector3& tan_point,
                                      const Vector3& normal) {
  Vector3 T = theta + radius * normal - tan_point;
  double K = (tan_point - theta).dot(normal) - radius;
  return alpha1 * T.dot(T) + alpha2 * std::pow(K, 2);
}
// sphere (theta, radiu)
double TangentPlane::update_energy_value(const Vector3& theta,
                                         const double& radius,
                                         const double alpha1,
                                         const double alpha2) {
  sanity_check();
  // Vector3 T = theta + radius * normal - points[0];
  // double K = (points[0] - theta).dot(normal) - radius;
  // energy = alpha1 * T.dot(T) + alpha2 * std::pow(K, 2);
  energy = get_energy_value(theta, radius, alpha1, alpha2, points[0], normal);
  energy_over_sq_radius = energy / std::pow(radius, 2);
  return energy;
}

void TangentPlane::copy_points_from(const TangentPlane& b) {
  std::copy(b.points.begin(), b.points.end(), std::back_inserter(this->points));
}

void TangentPlane::sanity_check() const {
  if (is_deleted) {
    logger().error("[TangentPlane] normal ({},{},{}) is_deleted {}", normal[0],
                   normal[1], normal[2], is_deleted);
    log_and_throw("ERROR");
  }
  if (points.size() != 1) {
    logger().error("[TangentPlane] normal ({},{},{}) points size {} > 1",
                   normal[0], normal[1], normal[2], points.size());
    log_and_throw("ERROR");
  }
}

void TangentPlane::update_extrude_ratio(const Vector3& center,
                                        const double& sq_radius) {
  double radius = std::sqrt(sq_radius);
  // double dist = std::abs(radius - (this->points[0] -
  // center).dot(this->normal));
  double d = get_dist_to_a_plane(this->points[0], this->normal, center);
  double dist = radius - d;
  this->extrude_ratio = std::abs(dist) / radius;
  this->extrude_ratio_signed = dist / radius;
}

/////////////////////////////////////////////////////////////////////////////////////
// ConnectedComponent
/////////////////////////////////////////////////////////////////////////////////////
void ConnectedComponent::clear() {
  adj_sphere_tags.clear();
  rpd_boundary_vertices.clear();
  rpd_facet_indices.clear();
  rpd_ref_fs_indices.clear();
}

void ConnectedComponent::clear_set_storage() {
  adj_sphere_tags.clear();
  rpd_boundary_vertices.clear();
  rpd_facet_indices.clear();
  rpd_ref_fs_indices.clear();
}

void ConnectedComponent::print_info() const {
  logger().debug("----- ConnectedComponent");
  logger().debug("is_deleted: {}, is_on_cc_line: {}", is_deleted,
                 is_on_cc_line);
  logger().debug(" rpd_facet_indices: {}, rpd_ref_fs_indices: {}",
                 rpd_facet_indices, rpd_ref_fs_indices);
  logger().debug("adj_sphere_tags: {}, rpd_boundary_vertices: {}",
                 adj_sphere_tags, rpd_boundary_vertices);
}

bool ConnectedComponent::is_adj_to(const int adj_tag) const {
  if (adj_sphere_tags.find(adj_tag) != adj_sphere_tags.end()) {
    return true;
  }
  return false;
}

/////////////////////////////////////////////////////////////////////////////////////
// MVertex
/////////////////////////////////////////////////////////////////////////////////////
MVertex::MVertex(int t, double v0, double v1, double v2, double sq_r,
                 SphereType _st) {
  tag = t;
  type = _st;
  pos = Vector3(v0, v1, v2);
  sq_radius = sq_r;
  dilate_sphere_radius();
  dt_seed_ids.resize(4);
  dt_seed_pos.resize(4);
  dt_seed_normals.resize(4);
  valid_idx = -1;
  matfp::color_addon_helper(color_rgb[0], color_rgb[1], color_rgb[2]);
}

MVertex::MVertex(int t, double v0, double v1, double sq_r) {
  tag = t;
  sq_radius = sq_r;
  pos_2d = Vector2(v0, v1);
}

// manually inserted new sphere for concave lines
MVertex::MVertex(const int t, const Vector3& pin, const Vector3& pin_normal) {
  tag = t;
  type = SphereType::T_2_c;
  ss.p = pin;
  ss.pnormal = pin_normal;
  double radius = INIT_RADIUS;
  pos = pin - pin_normal * radius;
  sq_radius = std::pow(radius, 2);
  is_outside = false;
}

void MVertex::print_info() const {
  logger().info("------ MVertex Info ------");
  logger().info(
      "tag: {}, valid_idx {}, is_deleted: {}, is_outside: {}, is_on_bbox: {}, "
      "is_a_feature_sphere: {}, is_on_corner: {}, min_corner_len: {}, "
      "is_a_concave_sphere {}, "
      "is_on_s_edge: {}, is_added_for_corner: {}, is_sphere_updated: {}",
      tag, valid_idx, is_deleted, is_outside, is_on_bbox, is_a_feature_sphere(),
      is_on_corner, min_corner_len, is_a_concave_sphere(), is_on_s_edge,
      is_added_for_corner, is_sphere_updated);
  logger().info(
      "pos: ({},{},{}) sq_radius: {}, radius: {}, sq_radius_dilated: {}, "
      "radius_dilated: {}",
      pos[0], pos[1], pos[2], sq_radius, std::sqrt(sq_radius),
      sq_radius_dilated, std::sqrt(sq_radius_dilated));
  logger().info("type: {}, added_for_two_spheres: {}, dt_seed_ids: {}", type,
                added_for_two_spheres, dt_seed_ids);
  logger().info("#cc_parts: {}, #tan_planes: {}, #tan_cc_lines: {}",
                cc_parts.size(), tan_planes.size(), tan_cc_lines.size());

  logger().info("ss_params: p: ({},{},{}), q: ({},{},{}), q_fidx: {}", ss.p[0],
                ss.p[1], ss.p[2], ss.q[0], ss.q[1], ss.q[2], ss.q_fid);

  if (is_a_feature_sphere()) logger().info("se_adj_se: {}", se_adj_se);
  // logger().debug("step_2_sq_energy: {}", step_2_sq_energy);
  logger().debug("newly_added_feature_mat_tags: {}",
                 newly_added_feature_mat_tags);
}

void MVertex::print_tan_planes() const {
  logger().debug("------- MVertex {} Tangent Planes: {}", tag,
                 tan_planes.size());
  for (const auto& tan_pl : tan_planes) {
    tan_pl.print_info();
  }
}

void MVertex::print_cc_parts() const {
  logger().debug("------- MVertex {} Connected Components: {}", tag,
                 cc_parts.size());
  for (const auto& tan_part : cc_parts) {
    tan_part.print_info();
  }
}

void MVertex::print_cc_lines() const {
  logger().debug("------- MVertex {} Concave Lines: {}", tag,
                 tan_cc_lines.size());
  for (const auto& cc_line : tan_cc_lines) {
    cc_line.print_info();
  }
}

/**
 * Compares two 4D points with respect to the lexicographic order.
 * s1, s2 are the coordinates of the two 4D points.
 * true if s1 is strictly before s2 in the lexicographic order.
 * false otherwise.
 */
// SCALAR_ZERO seems to small, use SCALAR_ZERO_N3 instead
bool MVertex::operator<(const MVertex& m2) const {
  if (std::abs(pos[0] - m2.pos[0]) > SCALAR_ZERO_N3) {
    return pos[0] < m2.pos[0];
  }
  if (std::abs(pos[1] - m2.pos[1]) > SCALAR_ZERO_N3) {
    return pos[1] < m2.pos[1];
  }
  if (std::abs(pos[2] - m2.pos[2]) > SCALAR_ZERO_N3) {
    return pos[2] < m2.pos[2];
  }
  if (std::abs(sq_radius - m2.sq_radius) > SCALAR_ZERO_N3) {
    return sq_radius < m2.sq_radius;
  }

  /*
   * Read this if you see heap overflow!!
   * Compare requirement (https://en.cppreference.com/w/cpp/named_req/Compare)
   * For all a, comp(a,a)==false
   * If comp(a,b)==true then comp(b,a)==false
   * otherwise the std::sort would return heap-buffer-overflow
   */

  if (is_deleted && !m2.is_deleted) {
    return true;  // // cur < m2, not deleted sphere larger
  }

  // check internal features sizes
  // keep the one with more internal feature neighbors
  if (intf_adj_intf.size() != m2.intf_adj_intf.size()) {
    return intf_adj_intf.size() < m2.intf_adj_intf.size();
  }
  // keep feature spheres large
  if (m2.type == SphereType::T_1_N) {
    if (type != SphereType::T_1_N) {
      return true;  // cur < m2
    }
  }
  if (type == SphereType::T_1_N) {
    if (m2.type != SphereType::T_1_N) {
      return false;  // m2 < cur
    }
  }

  // keep sphere crossing concave lines large
  if (m2.type == SphereType::T_2_c) {
    if (type != SphereType::T_2_c) return true;  // cur < m2
  }

  // also add type info, we want to keep the one with higher type
  if (type != SphereType::T_UNK && m2.type != SphereType::T_UNK) {
    return type < m2.type;
  }
  if (type > 0 && type > 0) {
    return type < m2.type;
  }

  // when a == b, then always return false to avoid heap overflow
  return false;
}

// for removing medial spheres that are too close
// using absolute value is not reliable for all models
// let's use ratio to measure how deep two spheres intersect
bool MVertex::operator==(const MVertex& m2) const {
  // our mesh has been normalized through
  // MeshIO::normalize_mesh()
  double threshold = 0.01;
  double diff = (pos - m2.pos).norm();
  if (diff > threshold) return false;
  diff = std::abs(std::sqrt(sq_radius) - std::sqrt(m2.sq_radius));
  if (diff > threshold) return false;
  return true;
}

// for sanity check
bool MVertex::operator!=(const MVertex& m2) const {
  if (tag != m2.tag) return true;
  double diff = (pos - m2.pos).norm();
  if (diff > SCALAR_ZERO_N3) return true;
  diff = std::abs(sq_radius - m2.sq_radius);
  if (diff > SCALAR_ZERO_N3) return true;
  return false;
}

bool MVertex::ensure_valid() {
  if (!pos.allFinite()) {
    logger().error(
        "[WRONG] mat {} has pos not all finite, delete and reset to 0!", tag);
    pos = Vector3::Zero();
    is_deleted = true;
    return false;
  }
  return true;
}

int MVertex::get_num_cc_parts_exclude_cc_lines() const {
  int count = 0;
  for (const auto& one_cc_part : cc_parts) {
    if (one_cc_part.is_on_cc_line) continue;
    count += 1;
  }
  return count;
}

bool MVertex::is_a_feature_sphere() const {
  if (type == SphereType::T_1_N) return true;
  if (is_on_s_edge || is_on_cc_edge) return true;
  if (is_on_corner) return true;
  return false;
}

bool MVertex::is_a_concave_sphere() const {
  if (type == SphereType::T_2_c || type == SphereType::T_N_c) {
    return true;
  }
  return false;
}

void MVertex::get_sphere_all_tangent_normals_includes_cc_lines(
    std::vector<Vector3>& normals) const {
  normals.clear();
  if (is_a_concave_sphere() || tan_cc_lines.size() > 0) {
    for (const auto& cc_line : tan_cc_lines) {
      normals.push_back(cc_line.adj_ref_normals[0]);
      normals.push_back(cc_line.adj_ref_normals[1]);
    }
  }
  for (const auto& tan_pl : tan_planes) {
    normals.push_back(tan_pl.normal);
  }
}

// tangent pair: (point, normal)
void MVertex::get_sphere_all_tangent_pairs_includes_cc_lines(
    std::vector<std::array<Vector3, 2>>& tan_pairs) const {
  tan_pairs.clear();
  if (is_a_concave_sphere() || tan_cc_lines.size() > 0) {
    for (const auto& cc_line : tan_cc_lines) {
      tan_pairs.push_back({{cc_line.tan_point, cc_line.adj_ref_normals[0]}});
      tan_pairs.push_back({{cc_line.tan_point, cc_line.adj_ref_normals[1]}});
    }
  }
  for (const auto& tan_pl : tan_planes) {
    tan_pairs.push_back({{tan_pl.points[0], tan_pl.normal}});
  }
}

// including spheres on sharp edges and corners
void MVertex::add_tan_planes_for_feature_sphere(
    const std::vector<Vector3>& normals,
    const std::vector<int>& ref_fids /* matching GEO::Mesh*/) {
  for (int i = 0; i < normals.size(); i++) {
    const Vector3 n = normals[i];
    const int fid = ref_fids[i];
    if (is_normal_covered_by_tan_pls(n)) continue;
    TangentPlane tan_pl(n);
    tan_pl.points.push_back(pos);
    tan_pl.fid = fid;
    tan_planes.push_back(tan_pl);
  }
}

bool MVertex::is_normal_covered_by_tan_pls(
    const Vector3& n, const double esp_degree_given) const {
  for (const auto& tan_pl : tan_planes) {
    if (tan_pl.is_same_normal(n, esp_degree_given)) {
      return true;
    }
  }
  return false;
}

bool MVertex::is_normal_covered_by_cc_lines(
    const Vector3& n, const double esp_degree_given) const {
  for (const auto& one_cc_line : tan_cc_lines) {
    if (one_cc_line.is_normal_covered_by_adj_fs(n, esp_degree_given)) {
      return true;
    }
  }
  return false;
}

bool MVertex::is_normal_covered_then_merge(const TangentPlane& tan_pl_new,
                                           const double esp_degree_given) {
  for (auto& tan_pl : tan_planes) {
    if (tan_pl.is_same_normal(tan_pl_new.normal, esp_degree_given)) {
      tan_pl.copy_points_from(tan_pl_new);
      return true;
    }
  }
  return false;
}

bool MVertex::is_normal_covered_includes_cc_lines(
    const Vector3& n, const double esp_degree_given) const {
  bool is_covered = is_normal_covered_by_tan_pls(n, esp_degree_given);
  for (const auto& one_cc_line : tan_cc_lines) {
    if (one_cc_line.is_normal_covered_by_adj_fs(n, esp_degree_given)) {
      is_covered = true;
      break;
    }
  }
  return is_covered;
}

void MVertex::update_tan_planes_points() {
  for (auto& tan_pl : tan_planes) {
    tan_pl.update_points_with_centroid();
  }
}

// mainly for updating fids
void MVertex::update_tan_planes_by_sf_mesh(const GEO::Mesh& sf_mesh,
                                           const AABBWrapper& aabb_wrapper) {
  // do not update for se_spheres
  if (is_a_feature_sphere()) {
    // check not fid == -1
    for (auto& tan_pl : tan_planes) {
      if (tan_pl.fid == -1) {
        logger().debug(
            "mat_p {} is a feature sphere, but has tangent plane with fid == "
            "-1",
            tag);
        log_and_throw("ERROR");
      }
    }
    return;
  }
  for (auto& tan_pl : tan_planes) {
    tan_pl.update_by_sf_mesh(sf_mesh, aabb_wrapper);
  }
}

bool MVertex::update_sphere_ss_param(bool is_debug) {
  if (is_outside || is_a_feature_sphere()) return false;

  bool is_param_saved = false;
  for (int i = 0; i < tan_planes.size(); i++) {
    const auto& tan_pl = tan_planes[i];
    if (tan_pl.is_deleted) continue;
    ss.p = tan_planes[i].points[0];
    ss.pnormal = tan_planes[i].normal;
    is_param_saved = true;
    break;
  }
  return is_param_saved;
}

// Initialize sphere to be tangent to pin point p with normal pnormal
// calling after ss_params are set
bool MVertex::init_ss_sphere(bool is_debug) {
  if (!ss.p.allFinite() || !ss.pnormal.allFinite()) return false;

  // Using the original radius from DT
  Vector3 center = ss.p - ss.pnormal * std::sqrt(sq_radius);
  if (!center.allFinite()) return false;

  if (is_debug) {
    logger().debug(
        "[ss_param] mat_p {} initialize sphere center from ({},{},{}) to "
        "({},{},{})...",
        tag, pos[0], pos[1], pos[2], center[0], center[1], center[2]);
  }
  pos = center;
  return true;
}

// Updating tangnet planes of sphere using ss:
// p and pnormal, q and qnormal
void MVertex::update_tan_planes_by_ss(bool is_debug) {
  if (is_debug) logger().debug("Updating tangent planes for mat_p: {}", tag);
  tan_planes.clear();

  if (is_a_concave_sphere()) {
    // the pin p is on concave line
    // do not save as tangent plane!
    if (tan_cc_lines.empty()) {
      print_info();
      log_and_throw("Concave sphere must have concave line");
    }
    tan_cc_lines[0].tan_point = ss.p;
    tan_cc_lines[0].normal = ss.pnormal;
  } else {
    TangentPlane tan_pl_p;
    tan_pl_p.points.push_back(ss.p);
    tan_pl_p.normal = ss.pnormal;
    tan_planes.push_back(tan_pl_p);
  }

  TangentPlane tan_pl_q;
  tan_pl_q.points.push_back(ss.q);
  tan_pl_q.normal = ss.qnormal;
  tan_planes.push_back(tan_pl_q);
}

void MVertex::save_old_pos_radius(bool is_clear) {
  old_pos = pos;
  old_sq_radius = sq_radius;
  if (is_clear) {
    pos << 0.0, 0.0, 0.0;
    sq_radius = 0.;
  }
}

// signed!!!
// +: sphere extrudes
// -: sphere not touching
double MVertex::dist_to_plane_signed(const Vector3& pl_point,
                                     const Vector3& pl_normal) {
  double radius = std::sqrt(this->sq_radius);
  // double dist = std::abs(radius - (pl_point - this->pos).dot(pl_normal));
  double d = get_dist_to_a_plane(pl_point, pl_normal, pos);
  double dist = radius - d;
  return dist;
}

double MVertex::get_extrude_ratio_unsigned(const Vector3& pl_point,
                                           const Vector3& pl_normal) {
  double radius = std::sqrt(sq_radius);
  double dist = std::abs(dist_to_plane_signed(pl_point, pl_normal));
  return dist / radius;
}
double MVertex::get_extrude_ratio_signed(const Vector3& pl_point,
                                         const Vector3& pl_normal) {
  double radius = std::sqrt(sq_radius);
  double dist = dist_to_plane_signed(pl_point, pl_normal);
  return dist / radius;
}

double MVertex::max_dist_to_tan_points(
    const std::vector<TangentPlane>& given_tan_pls) {
  double max_dist = DBL_MIN;
  for (const auto& tan_pl : given_tan_pls) {
    double dist = sphere_dist_to_point(pos, std::sqrt(sq_radius),
                                       tan_pl.points[0], tan_pl.normal);
    if (dist > max_dist) {
      max_dist = dist;
    }
  }
  return max_dist;
}

// Dilate the sphere radius to make it more robust
void MVertex::dilate_sphere_radius() {
  if (this->is_a_feature_sphere()) {
    this->sq_radius_dilated = std::pow(0.01, 2);
    return;
  }
  double radius = std::sqrt(this->sq_radius);
  if (type == SphereType::T_N || type == SphereType::T_N_c ||
      is_added_for_corner) {
    // more important
    radius += radius * 0.06;
  } else {
    radius += radius * 0.05;
  }
  this->sq_radius_dilated = std::pow(radius, 2);
}

double MVertex::update_all_energy_values(const double alpha1,
                                         const double alpha2,
                                         const double alpha3) {
  double radius = std::sqrt(sq_radius);
  double sum = 0.;
  for (auto& tan_pl : tan_planes) {
    tan_pl.update_energy_value(pos, radius, alpha1, alpha2);
    sum += tan_pl.energy;
  }
  for (auto& one_cc_line : tan_cc_lines) {
    one_cc_line.update_energy_value(pos, radius, alpha3);
    sum += one_cc_line.energy;
  }
  return sum;
}

void MVertex::store_sq_energy(const AABBWrapper& aabb_wrapper) {
  double sq_dist = aabb_wrapper.get_sq_dist_to_sf(pos);
  double sq_energy = std::pow(std::sqrt(sq_dist) - std::sqrt(sq_radius), 2);
  // logger().debug("mat_p {} has sq_dist {}, sq_radius: {}, sq_energy: {}",
  // tag,
  //                sq_dist, sq_radius, sq_energy);
  step_2_sq_energy.push_back(sq_energy);
}

/////////////////////////////////////////////////////////////////////////////////////
// Corner Preservation
/////////////////////////////////////////////////////////////////////////////////////
void SHEET_SEG::print_info() const {
  logger().debug(
      "[Sheet Segment]: has {} tangent planes, {} tangent concave lines: {}",
      tan_pls.size(), tan_cc_line_ids.size(), tan_cc_line_ids);
  for (const auto& tan_pl : tan_pls) {
    tan_pl.print_info();
  }
}

void SHEET_NODE::print_info() const {
  logger().debug("-------- Sheet Node id {}:", id);
  logger().debug(
      "type: {}, related_mat_tag: {}, corner_tag: {}, sp_edge_id: {}", type,
      related_mat_tag, corner_tag, sp_edge_id);
  for (const auto& seg : sheet_segs) {
    seg.print_info();
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Export
/////////////////////////////////////////////////////////////////////////////////////
void export_medial_spheres(const std::string& maname,
                           const std::vector<MVertex>& all_medial_spheres) {
  std::string ma_name_full = maname + "_spheres_step1" + ".ma";
  logger().debug("saving medial spheres to .ma file: {}", ma_name_full);

  std::ofstream fout(ma_name_full);
  fout << all_medial_spheres.size() << std::endl;

  for (int i = 0; i < all_medial_spheres.size(); i++) {
    const auto& mat_p = all_medial_spheres[i];
    if (mat_p.is_outside) continue;
    // if (mat_p.valid_idx == -1)
    // 	continue;
    Vector3 pos = mat_p.pos;
    double radius = std::sqrt(mat_p.sq_radius);
    fout << "v " << pos[0] << " " << pos[1] << " " << pos[2] << " " << radius
         << " " << int(mat_p.type) << " " << mat_p.is_deleted_int << " "
         << mat_p.is_addition_int << " " << std::endl;
  }

  fout.close();
}

// format:
// s t1 t2 t3 ...
// step1 avg_e_type1 avg_e_type2 avg_e_type3 ...
// step2 avg_e_type1 avg_e_type2 avg_e_type3 ...
// step2 avg_e_type1 avg_e_type2 avg_e_type3 ...
void export_sphere_energy_spline(
    const std::string& maname, const std::vector<MVertex>& all_medial_spheres) {
  // std::string ma_name_full =
  //     "../out/opt/opt_" + maname + "_" + get_timestamp() + ".txt";
  std::string ma_name_full = "../out/opt/opt_" + maname + ".txt";
  logger().debug("start saving optimization .txt file: {}", ma_name_full);

  // step -> type -> sq_energy
  std::map<int, std::map<int, double>> step_2_type_2_sq_energies;
  std::map<int, int> type_2_cnt;  // type -> number of spheres
  for (int i = 0; i < all_medial_spheres.size(); i++) {
    const auto& mat_p = all_medial_spheres[i];
    const int type = mat_p.type;
    if (mat_p.is_outside || mat_p.is_deleted) continue;
    if (type == SphereType::T_UNK) continue;
    if (mat_p.is_a_feature_sphere()) continue;
    // merge for each step
    // logger().debug("mat_p {} with type {} has energies: {}", mat_p.tag,
    //                mat_p.type, mat_p.step_2_sq_energy);
    for (int step = 0; step < mat_p.step_2_sq_energy.size(); step++) {
      auto& type_2_sq_e = step_2_type_2_sq_energies[step];
      type_2_sq_e[mat_p.type] += mat_p.step_2_sq_energy[step];
    }
    auto& type_cnt = type_2_cnt[mat_p.type];
    type_cnt++;
  }
  // logger().debug("before average step_2_type_2_sq_energies: {}",
  //  step_2_type_2_sq_energies);

  // count average
  for (auto& pair : step_2_type_2_sq_energies) {
    int step = pair.first;
    auto& type_2_sq_energies = pair.second;
    for (auto& pair2 : type_2_sq_energies) {
      int type = pair2.first;
      auto& sq_energy = pair2.second;
      int cnt = type_2_cnt.at(type);
      sq_energy /= cnt;
    }
  }
  logger().debug("step_2_type_2_sq_energies: {}", step_2_type_2_sq_energies);

  std::vector<int> all_types = {{2, 3, 11, 12}};
  std::ofstream fout(ma_name_full);
  for (auto& pair : step_2_type_2_sq_energies) {
    int step = pair.first;
    fout << step;
    auto& type_2_sq_energies = pair.second;
    for (int type : all_types) {
      if (type_2_sq_energies.find(type) != type_2_sq_energies.end()) {
        const auto& sq_e = type_2_sq_energies.at(type);
        fout << " " << sq_e;
      } else {
        fout << " " << 0;
      }
    }
    fout << std::endl;
  }
  fout.close();
}

}  // namespace matfp