// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/ShrinkSpheres.h"

#include "matfp/InscribedSpheres.h"

namespace matfp {

bool is_sphere_need_to_be_updated_copy(MVertex& mat_p) {
  bool is_update = false;
  // check all tangent planes, if extrude too much, then delete
  const Vector3& center = mat_p.pos;
  const double radius = std::sqrt(mat_p.sq_radius);
  for (auto& tan_pl : mat_p.tan_planes) {
    tan_pl.update_extrude_ratio(mat_p.pos, mat_p.sq_radius);
    if (!is_update && tan_pl.extrude_ratio > extrude_threshold) {
      is_update = true;
      logger().debug(
          "[IS_UPDATE] mat_p {} has tan_pl.extrude_ratio {}, going to be "
          "updated",
          mat_p.tag, tan_pl.extrude_ratio);
      // do not break here
      // we want to update all extrude ratio
    }
  }
  return is_update;
}

// plane is defined by point p and normal pn
bool is_sphere_extrude_plane(const Vector3& center, const double& radius,
                             const Vector3& p, const Vector3& pn) {
  // Calculate new q
  double d = get_dist_to_a_plane(p, pn, center);
  double diff_ratio = std::abs(radius - d) / radius;
  if (diff_ratio <= extrude_threshold_2) {
    return false;
  }
  return true;
}

void update_sphere_ss_param_using_seeds(
    const int nb_seeds, const std::set<int>& feature_points,
    std::vector<MVertex>& all_medial_spheres) {
  // for selecting pin point of ShrinkSpheres
  std::vector<bool> is_pts_selected(nb_seeds, false);

  // if feature points, then we set true for ShrinkSphere
  for (const auto& i : feature_points) {
    is_pts_selected[i] = false;
  }

  // set pin point for each inner non-feature medial sphere
  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_outside || mat_p.is_a_feature_sphere()) continue;

    // for ShrinkSpheres
    bool is_found_pin = false;
    for (int i = 0; i < 4; i++) {
      int pid = mat_p.dt_seed_ids[i];
      if (pid > nb_seeds)
        // could be 8 bbox points?
        continue;

      if (!is_pts_selected[pid]) {
        mat_p.ss.pidx = pid;
        mat_p.ss.p = mat_p.dt_seed_pos[i];
        mat_p.ss.pnormal = mat_p.dt_seed_normals[i];
        is_pts_selected[pid] = true;
        is_found_pin = true;
        break;
      }
    }
    if (!is_found_pin) {
      logger().error("sphere tag {} did not find pin point!!!", mat_p.tag);
    }
  }
}

double compute_radius(const Vector3& p, const Vector3& n, const Vector3& q,
                      bool is_debug = false) {
  // Compute radius of the ball that touches points p and q and whose center
  // falls on the normal n from p
  Scalar d = (p - q).norm();
  Scalar cos_theta = n.dot(p - q) / d;

  if (is_debug)
    logger().debug(
        "p: ({},{},{}), q: ({},{},{}), n: ({},{},{}), d: {}, cos_theta: {}",
        p[0], p[1], p[2], q[0], q[1], q[2], n[0], n[1], n[2], d, cos_theta);

  return Scalar(d / (2 * cos_theta));
}

// Calculate a medial ball for a given oriented point using the shrinking ball
// algorithm, see https://3d.bk.tudelft.nl/rypeters/pdfs/16candg.pdf section 3.2
// for details
bool shrink_sphere(const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                   MVertex& mat_p, bool is_debug = false) {
  const unsigned int iteration_limit = 30;
  const Vector3& p = mat_p.ss.p;
  const Vector3& n = mat_p.ss.pnormal;
  // copy
  Vector3 c_next = mat_p.pos;
  double r_next = std::sqrt(mat_p.sq_radius);
  double sq_dist = -1;
  int num_itr = 0;
  // mat_p.store_sq_energy(aabb_wrapper);  // init
  bool is_good = true;

  while (true) {
    int q_fid_prev_prev = mat_p.ss.q_fid_prev;
    int q_fid_prev = mat_p.ss.q_fid;
    Vector3 q;
    int q_fid = aabb_wrapper.get_nearest_point_on_sf(c_next, q, sq_dist);
    const Vector3 qnormal = get_mesh_facet_normal(sf_mesh, q_fid);
    if (is_debug) {
      logger().debug("[Shrink] mat_p {} has q_fid: {}->{}->{}", mat_p.tag,
                     q_fid_prev_prev, q_fid_prev, q_fid);
    }

    // This should handle all (special) cases where we want to break the loop
    // - normal case when ball no longer shrinks
    // - the case where q==p
    // - any duplicate point cases
    // Detail see:
    // https://github.com/tudelft3d/masbcpp/blob/master/src/compute_ma_processing.cpp#L88
    double dist_pq = get_distance_between_two_vectors(p, q);
    double sq_radius_eps =
        (r_next - SCALAR_ZERO_N3) * (r_next - SCALAR_ZERO_N3);
    if (sq_radius_eps < SCALAR_ZERO_N4) {
      if (is_debug)
        logger().error("[Shrink] mat_p {} shrink too small, delete", mat_p.tag);
      is_good = false;
      break;
    }
    if (sq_dist >= sq_radius_eps || dist_pq < SCALAR_ZERO_N3) {
      if (is_debug)
        logger().debug("[Shrink] break, mat_p {} found its breaking condition",
                       mat_p.tag);
      break;
    }

    bool is_extrude = is_sphere_extrude_plane(c_next, r_next, q, qnormal);
    if (!is_extrude) {
      if (is_debug)
        logger().debug("[Shrink] break, mat_p {} not extrude from q",
                       mat_p.tag);
      break;
    }

    // Compute next ball center
    r_next = compute_radius(p, n, q, is_debug);
    c_next = p - n * r_next;

    if (!c_next.allFinite()) {
      if (is_debug)
        logger().debug(
            "[Shrink] break, mat_p {} has q_fid: {}, c_next ({},{},{}) "
            "with "
            "r_next: {}",
            mat_p.tag, q_fid, c_next[0], c_next[1], c_next[2], r_next);
      break;
    }

    // Stop iteration if this looks like an infinite loop:
    if (num_itr > iteration_limit) {
      if (is_debug)
        logger().debug("[Shrink] break, mat_p {} reaches iteration limits",
                       mat_p.tag);
      break;
    }

    if (is_debug) {
      logger().debug(
          "[Shrink] mat_p {} updating, pos ({},{},{}) -> ({},{},{}), radius "
          "{}->{}. q_fid: {}->{}->{}",
          mat_p.tag, mat_p.pos[0], mat_p.pos[1], mat_p.pos[2], c_next[0],
          c_next[1], c_next[2], std::sqrt(mat_p.sq_radius), r_next,
          q_fid_prev_prev, q_fid_prev, q_fid);
    }

    mat_p.pos = c_next;
    mat_p.sq_radius = std::pow(r_next, 2);
    mat_p.ss.q = q;
    mat_p.ss.qnormal = qnormal;
    mat_p.ss.q_fid_prev = mat_p.ss.q_fid;
    mat_p.ss.q_fid = q_fid;
    num_itr++;
    // mat_p.store_sq_energy(aabb_wrapper);
  }  // while true
  if (is_debug) {
    logger().debug("mat_p: {} has num_itr: {}", mat_p.tag, num_itr);
  }
  if (is_good == false || num_itr == 0) {
    return false;
  }
  return true;
}

// a wrapper for calling shrink_sphere()
// could be SphereType::T_2 or SphereType::T_2_c
bool shrink_sphere_wrapper(const GEO::Mesh& sf_mesh,
                           const AABBWrapper& aabb_wrapper, MVertex& mat_p,
                           bool is_setup_ss_param, bool is_check_cc,
                           bool is_debug) {
  bool is_success = false;
  // For sphere that not setup ss_param
  // eg. spheres inserted for concave lines
  if (is_setup_ss_param) {
    is_success = mat_p.update_sphere_ss_param();
    if (!is_success) {
      // logger().error("ERROR: mat {} not update ss parameters!!!", mat_p.tag);
      return false;
    }
    is_success = mat_p.init_ss_sphere(is_debug);
    if (!is_success) {
      // logger().error("ERROR: mat {} not able to initialize!!!", mat_p.tag);
      return false;
    }
  }

  is_success = shrink_sphere(sf_mesh, aabb_wrapper, mat_p, is_debug);
  if (!is_success) {
    if (is_debug)
      logger().error("ERROR: mat {} not shrink properly, delete", mat_p.tag);
    return false;
  }
  // update mat_p.tan_planes by p and q from shrinking sphere
  mat_p.update_tan_planes_by_ss();
  if (mat_p.type != SphereType::T_2_c) {
    mat_p.type = SphereType::T_2;
  }
  mat_p.dilate_sphere_radius();

  // all in false for now
  if (is_check_cc) {
    // post checking T_2 spheres
    // if the p/q point is too close to concave lines, then delete
    // this is because we have inserted enough spheres around concave lines
    // by setting pin point p on concave lines and shrink
    double p_sq_dist_2cc = aabb_wrapper.get_sq_dist_to_cc(mat_p.ss.p);
    double q_sq_dist_2cc = aabb_wrapper.get_sq_dist_to_cc(mat_p.ss.q);
    if (mat_p.type != SphereType::T_2_c &&
        (p_sq_dist_2cc <= SQ_DIST_TO_CC || q_sq_dist_2cc <= SQ_DIST_TO_CC)) {
      // logger().error("ERROR: mat {} shrink too close to concave lines,
      // delete",
      //                mat_p.tag);
      return false;
    }
  }

  return true;
}

// For all T_2 spheres
// will remove T_2 spheres that close to concave lines
void shrink_spheres(const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                    const std::vector<TangentConcaveLine>& tan_cc_lines,
                    std::vector<MVertex>& all_medial_spheres, bool is_debug) {
  for (MVertex& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) continue;
    if (mat_p.is_on_s_edge || mat_p.is_on_cc_edge) continue;
    if (mat_p.is_on_corner || mat_p.is_added_for_corner) continue;

    // only update newly inserted spheres for concave lines
    // if (mat_p.type != SphereType::T_2_c) continue;

    // update T_2_c and T_2 spheres
    if (mat_p.type != SphereType::T_2_c && mat_p.tan_planes.size() != 2)
      continue;
    mat_p.is_sphere_updated = true;  // for add_or_delete_for_se()

    if (is_debug) {
      logger().debug("[SHRINK] start shrinking sphere mat_p {} ...", mat_p.tag);
    }

    // For spheres from DT
    bool is_setup_ss_param = true;
    if (mat_p.type == SphereType::T_2_c) {
      is_setup_ss_param = false;  // ss_param already setup
    }
    bool is_success =
        shrink_sphere_wrapper(sf_mesh, aabb_wrapper, mat_p, is_setup_ss_param,
                              true /*is_check_cc*/, is_debug);
    if (!is_success) {
      mat_p.is_deleted = true;
      continue;
    }

    // Make sphere of type SphereType::T_N_c
    // check if q is also around concave line
    // if so, then remove from tangent plane and save as concave line
    // (sphere tangent to 2 concave lines)
    if (mat_p.type == SphereType::T_2_c) {
      double q_sq_dist;
      Vector3 q_tan_point = mat_p.ss.q;
      int eid =
          aabb_wrapper.project_to_cc_get_nearest_face(q_tan_point, q_sq_dist);
      double threshold_over_sq_radius = q_sq_dist / mat_p.sq_radius;
      if (threshold_over_sq_radius <= SQ_DIST_TO_CC) {
        // only 1 tangent plane
        if (mat_p.tan_planes.size() != 1) {
          logger().error("[SHRINK] mat_p {} has {} tangent planes != 1",
                         mat_p.tag, mat_p.tan_planes.size());
          log_and_throw("ERROR");
        }
        TangentConcaveLine new_cc_line = tan_cc_lines[eid];
        new_cc_line.tan_point = q_tan_point;
        new_cc_line.normal = (q_tan_point - mat_p.pos).normalized();
        mat_p.tan_cc_lines.push_back(new_cc_line);
        mat_p.tan_planes.clear();
        // mat_p.type = SphereType::T_N_c;
        mat_p.type = SphereType::T_2_c;
        if (is_debug)
          logger().debug(
              "create T_N_c sphere {} with cc_line id: {}, "
              "threshold_over_sq_radius: {}",
              mat_p.tag, eid, threshold_over_sq_radius);
      }
    }  // SphereType::T_N_c

  }  // for all_medial_spheres

  //////////
  // remove duplicated spheres
  std::vector<MVertex> sorted_partial_medial_spheres;
  remove_duplicated_medial_spheres(sorted_partial_medial_spheres,
                                   all_medial_spheres, true /*is_override*/);
}

// create new spheres (given num_new_spheres)
// each with a given pin and pin_normal
//
// NOTE: num_new_spheres >= 2, including 2 boundary normals of concave line
int create_new_spheres_for_cc_line(const Vector3& pin,
                                   const TangentConcaveLine& one_cc_line,
                                   const double cc_normal_eps,
                                   std::vector<MVertex>& all_medial_spheres) {
  const std::array<Vector3, 2>& adj_ref_normals = one_cc_line.adj_ref_normals;
  const double angle = angle_between_two_vectors_in_degrees(adj_ref_normals[0],
                                                            adj_ref_normals[1]);
  // 2 boundary normals
  int num_new_spheres = 2;
  // new normals in between
  num_new_spheres += std::ceil(angle / cc_normal_eps) - 1;
  // logger().debug(
  //     "[CC Sphere] calling create_new_spheres_for_cc_line with "
  //     "num_new_spheres: {} =  "
  //     "angle: {} / eps: "
  //     "{} + 2",
  //     num_new_spheres, angle, cc_normal_eps);
  std::vector<Vector3> new_normals;  // size will be num_new_spheres
  std::vector<int> new_sphere_tags;
  sample_k_vectors_given_two_vectors(adj_ref_normals[0], adj_ref_normals[1],
                                     num_new_spheres, new_normals);
  for (int i = 0; i < num_new_spheres; i++) {
    MVertex new_sphere(all_medial_spheres.size(), pin, new_normals[i]);
    new_sphere.tan_cc_lines.push_back(one_cc_line);
    all_medial_spheres.push_back(new_sphere);
    new_sphere_tags.push_back(new_sphere.tag);
  }
  // logger().debug("[CC Sphere] cc_line {} creating new cc sphere: {}",
  //                one_cc_line.id, new_sphere_tags);
  return num_new_spheres;
}

// Insert new spheres that pin at concave edges.
// We specify pin points on concave edge with length cc_len_eps (user defined),
// and for each pin point, we specify the number of newly sampled normals with
// cc_normal_eps (user defined)
void insert_spheres_for_concave_lines(
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    std::vector<MVertex>& all_medial_spheres,
    double cc_len_eps, /*=length, scaled in [0,10]*/
    double cc_normal_eps /*=degree, scaled in [0, 360]*/) {
  // store if two end points of given one_cc_line is visited as pin
  // logger().debug(
  //     "[CC Sphere] tan_cc_lines: {}, cc_len_eps: {}, cc_normal_eps: {} ",
  //     tan_cc_lines.size(), cc_len_eps, cc_normal_eps);
  int num_new_spheres = 0;
  std::set<int> is_end_point_visited;
  for (const auto& one_cc_line : tan_cc_lines) {
    // for two end points
    for (int i = 0; i < 2; i++) {
      const int end_vs_id =
          one_cc_line.ref_vs_seed_ids[i];  // only fetch the first vs_pair
      const Vector3& pin = one_cc_line.ref_vs_pos[i];
      // logger().debug(
      //     "[CC Sphere] checking one_cc_line {} with end_vs_id: {}, pin: "
      //     "({},{},{})",
      //     one_cc_line.id, end_vs_id, pin[0], pin[1], pin[2]);
      // visted as pin already, skip
      if (is_end_point_visited.find(end_vs_id) != is_end_point_visited.end())
        continue;
      is_end_point_visited.insert(end_vs_id);
      num_new_spheres += create_new_spheres_for_cc_line(
          pin, one_cc_line, cc_normal_eps, all_medial_spheres);
    }

    // add new pin points given cc_len_eps
    const Vector3& start_pos = one_cc_line.ref_vs_pos[0];
    const Vector3& end_pos = one_cc_line.ref_vs_pos[1];
    Vector3 dir = get_direction(start_pos, end_pos);
    int step = std::floor(get_distance_between_two_vectors(end_pos, start_pos) /
                          cc_len_eps);
    // logger().debug(
    //     "[CC Sphere] cc_line {} has start_pos ({},{},{}), end_pos:
    //     ({},{},{}), " "dir: ({},{},{})", one_cc_line.id, start_pos[0],
    //     start_pos[1], start_pos[2], end_pos[0], end_pos[1], end_pos[2],
    //     dir[0], dir[1], dir[2]);
    for (int j = 1; j < step; j++) {
      Vector3 new_pin = start_pos + j * cc_len_eps * dir;
      // logger().debug(
      //     "[CC Sphere] cc_line {} created new_pin: ({},{},{}) using step
      //     {}/{}", one_cc_line.id, new_pin[0], new_pin[1], new_pin[2], j,
      //     step);
      num_new_spheres += create_new_spheres_for_cc_line(
          new_pin, one_cc_line, cc_normal_eps, all_medial_spheres);
    }
  }

  logger().debug("[CC Sphere] Created {} CC spheres", num_new_spheres);
}

}  // namespace matfp