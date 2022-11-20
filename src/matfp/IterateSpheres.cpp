// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/IterateSpheres.h"

#include <assert.h>

#include <Eigen/Dense>

#include "matfp/Common.h"
#include "matfp/InscribedSpheres.h"

namespace matfp {
/////////////////////////////////////////////////////////////////////////////////////
// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////
// This energy function has 3 parts with weights:
// 1. distance to tangent point -> alpha 1 (pN, nN)
// 2. distance to tangent plane -> alpha 2 (pN, nN)
// 3. distance to concave line -> alpha 3 (c_pN, c_nN)
bool update_sphere_given_plN_opt(
    const double alpha1, const double alpha2, const std::vector<Vector3>& pN,
    const std::vector<Vector3>& nN, const double alpha3,
    const std::vector<Vector3>& c_pN, const std::vector<Vector3>& c_nN,
    Vector3& new_center, double& new_sq_radius, bool is_debug) {
  if (pN.size() != nN.size() || c_pN.size() != c_nN.size()) {
    logger().error(
        "[TAN_PLN] unmatched plane parameters, pN: {}, nN: {}, c_pN: {}, c_nN: "
        "{}",
        pN.size(), nN.size(), c_pN.size(), c_nN.size());
    return false;
  }

  if (is_debug) {
    logger().debug(
        "[TAN_PLN] with N {}, cN {}, alpha1: {}, alpha2: {}, alpha3: {}",
        pN.size(), c_pN.size(), alpha1, alpha2, alpha3);
  }
  Matrix3 A = Matrix3::Zero();
  Vector3 B = Vector3::Zero();
  Vector3 C = Vector3::Zero();
  Vector3 E = Vector3::Zero();
  double D = 0., F = 0.;
  // if any tangent plane
  for (int i = 0; i < pN.size(); i++) {
    const auto& p = pN[i];
    const auto& n = nN[i];
    Matrix3 common =
        alpha1 * Matrix3::Identity() + alpha2 * (n * n.transpose());
    A += common;
    B += (alpha1 + alpha2) * n;
    C += (alpha1 + alpha2) * n;
    D += alpha1 + alpha2;
    E += common * p;
    F += (alpha1 + alpha2) * n.dot(p);
  }
  // if any concave line
  for (int i = 0; i < c_pN.size(); i++) {
    const auto& p = c_pN[i];
    const auto& n = c_nN[i];
    A += alpha3 * Matrix3::Identity();
    B += alpha3 * n;
    C += alpha3 * n;
    D += alpha3;
    E += alpha3 * p;
    F += alpha3 * n.dot(p);
  }

  Matrix4 N;
  N << A(0, 0), A(0, 1), A(0, 2), B[0], A(1, 0), A(1, 1), A(1, 2), B[1],
      A(2, 0), A(2, 1), A(2, 2), B[2], C[0], C[1], C[2], D;
  Vector4 b;
  b << E[0], E[1], E[2], F;

  Matrix4 N_inv;
  bool invertible;
  N.computeInverseWithCheck(N_inv, invertible);
  if (!invertible) {
    // logger().error("[TAN_PLN] Matrix 4x4 not invertible");
    // logger().debug(
    //     "[TAN_PLN] Matrix 4x4: \n"
    //     "({},{},{},{})\n({},{},{},{})\n({},{},{},{})\n({},{},{},{})\n",
    //     N(0, 0), N(0, 1), N(0, 2), N(0, 3), N(1, 0), N(1, 1), N(1, 2), N(1,
    //     3), N(2, 0), N(2, 1), N(2, 2), N(2, 3), N(3, 0), N(3, 1), N(3, 2),
    //     N(3, 3));
    return false;
  }

  // update mat_p
  Vector4 c_r = N_inv * b;
  if (is_debug) {
    logger().debug("[TAN_PLN] mat new sphere: center ({},{},{}), radius: {}",
                   c_r[0], c_r[1], c_r[2], c_r[3]);
  }
  if (c_r[3] < SCALAR_ZERO_N3) return false;
  new_center = Vector3(c_r[0], c_r[1], c_r[2]);
  new_sq_radius = std::pow(c_r[3], 2);
  return true;
}

// Given fid_given, we fetch the k_ring neighboring faces
// not crossing 1. sharp edges or 2. concave edges
// info stored in ref_fs_pairs_not_cross
void get_k_ring_neighbors_no_cross(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const int fid_given, const int k, std::set<int>& k_ring_fids,
    bool is_debug) {
  const GEO::Attribute<bool> v_is_on_se(sf_mesh.vertices.attributes(),
                                        "is_on_se");
  auto is_skip_se_neighbor = [&](const int f, const int nf) {
    std::array<int, 2> ref_fs_pair = {{f, nf}};
    std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
    if (ref_fs_pairs_not_cross.find(ref_fs_pair) !=
        ref_fs_pairs_not_cross.end()) {
      if (is_debug)
        logger().debug("[K_RING] face {} skip nf {} since sharing a sharp edge",
                       f, nf);
      return true;  // skip checking its neighbor
    }
    return false;
  };

  if (fid_given < 0) {
    logger().debug("[K_RING] error fid_given: {}", fid_given);
    log_and_throw("ERROR");
  }
  if (is_debug)
    logger().debug("calling get_k_ring_neighbor_no_se for fid_given: {}",
                   fid_given);

  k_ring_fids.clear();
  k_ring_fids.insert(fid_given);
  std::set<int> new_added_fids, new_added_fids_copy;
  new_added_fids.insert(fid_given);

  int i = 0;
  while (i < k) {
    new_added_fids_copy = new_added_fids;
    new_added_fids.clear();
    for (const auto& fid : new_added_fids_copy) {
      for (GEO::index_t le = 0; le < sf_mesh.facets.nb_vertices(fid); le++) {
        GEO::index_t nfid = sf_mesh.facets.adjacent(fid, le);
        if (nfid == GEO::NO_FACET) continue;
        if (is_skip_se_neighbor(fid, nfid)) continue;
        if (is_debug)
          logger().debug("fid {} has neighbor face nfid {}", fid, nfid);
        k_ring_fids.insert(nfid);
        new_added_fids.insert(nfid);
      }  // for facets.nb_vertices
    }    // for k_ring_fids
    i++;
  }  // while (i < k)
  if (is_debug) {
    logger().debug("[K_RING] fid {} found k_ring {} neighbors: {}", fid_given,
                   k, k_ring_fids);
  }
}

void get_faces_params(const GEO::Mesh& sf_mesh, std::set<int>& k_ring_fids,
                      std::map<int, std::array<Vector3, 4>>& faces_params) {
  // face_pts_normal stores a vector of triangle
  // each triangle we store fid -> three points (pA, pB, pC) and normal nN
  faces_params.clear();
  for (const int& fid : k_ring_fids) {
    std::array<Vector3, 4> params;
    for (int lv = 0; lv < 3; lv++) {
      params[lv] =
          to_eigen(sf_mesh.vertices.point(sf_mesh.facets.vertex(fid, lv)));
    }
    params[3] = get_mesh_facet_normal(sf_mesh, fid);
    faces_params[fid] = params;
  }
}

// Given sphere (theta, r), find the closest tangent point pN
// on triangle (pA, pB, pC) with normal nN
// by mininumizing energy function:
//
// This energy function has 2 parts with weights:
// 1. distance to tangent point -> alpha 1
// 2. distance to tangent plane -> alpha 2
double get_closest_tangent_point_opt(const double alpha1, const double alpha2,
                                     const Vector3& pA, const Vector3& pB,
                                     const Vector3& pC, const Vector3& nN,
                                     const Vector3& theta,
                                     const double& sq_radius, Vector3& pN,
                                     bool is_debug) {
  double radius = std::sqrt(sq_radius);
  Vector3 pAB = pB - pA;
  Vector3 pAC = pC - pA;
  Matrix3 N = nN * nN.transpose();  // symmetric matrix
  Vector3 M1 = theta + radius * nN;
  Vector3 M2 = N * theta + radius * nN;

  // NOTE:
  // here N is a symmetric matrix
  // so "pAB.transpose() * N * pAB" (not compiling) can be written
  // as "(N * pAB).dot(pAB)", otherwise its wrong
  Matrix2 K;
  K(0, 0) = alpha1 * (pAB.dot(pAB)) + alpha2 * (N * pAB).dot(pAB);
  K(0, 1) = alpha1 * (pAB.dot(pAC)) + alpha2 * (N * pAB).dot(pAC);
  K(1, 0) = alpha1 * (pAB.dot(pAC)) + alpha2 * (N * pAB).dot(pAC);
  K(1, 1) = alpha1 * (pAC.dot(pAC)) + alpha2 * (N * pAC).dot(pAC);
  Vector2 b;
  b[0] = (alpha1 * M1 + alpha2 * M2).dot(pAB) - alpha1 * pA.dot(pAB) -
         alpha2 * (N * pA).dot(pAB);
  b[1] = (alpha1 * M1 + alpha2 * M2).dot(pAC) - alpha1 * pA.dot(pAC) -
         alpha2 * (N * pA).dot(pAC);

  Matrix2 K_inv;
  bool invertible;
  K.computeInverseWithCheck(K_inv, invertible);
  if (!invertible) {
    if (is_debug) {
      logger().error("[CLOSEST_PNT] Matrix 2x2 not invertible");
      logger().debug("{}, {} \n {}, {}", K(0, 0), K(0, 1), K(1, 0), K(1, 1));
    }
    return DBL_MAX;
  }
  // update mat_p
  VectorXs t1_t2 = K_inv * b;
  double sum = t1_t2[0] + t1_t2[1];
  if (is_debug) {
    logger().debug("[CLOSEST_PNT] t1: {}, t2: {}, sum: {}", t1_t2[0], t1_t2[1],
                   sum);
  }
  if (sum <= 1 && sum >= 0 && t1_t2[0] <= 1 && t1_t2[0] >= 0 && t1_t2[1] <= 1 &&
      t1_t2[1] >= 0) {
    pN = pA + t1_t2[0] * pAB + t1_t2[1] * pAC;
    double E =
        TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pN, nN);
    if (is_debug)
      logger().debug("[CLOSEST_PNT] return pN ({},{},{}) with energy E {}",
                     pN[0], pN[1], pN[2], E);
    return E;
  }

  // closest point not in triangle, check 3 vertices
  double EA =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pA, nN);
  double EB =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pB, nN);
  double EC =
      TangentPlane::get_energy_value(theta, radius, alpha1, alpha2, pC, nN);
  if (EA <= EB && EA <= EC) {
    if (is_debug)
      logger().debug("[CLOSEST_PNT] return pA ({},{},{}) with energy EA {}",
                     pA[0], pA[1], pA[2], EA);
    pN = pA;
    return EA;
  }
  if (EB <= EA & EB <= EC) {
    if (is_debug)
      logger().debug("[CLOSEST_PNT] return pB ({},{},{}) with energy EB {}",
                     pB[0], pB[1], pB[2], EB);
    pN = pB;
    return EB;
  }

  if (is_debug)
    logger().debug("[CLOSEST_PNT] return pC ({},{},{}) with energy EC {}",
                   pC[0], pC[1], pC[2], EC);
  pN = pC;
  return EC;
}

void update_tangent_points_on_tan_pls(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross, MVertex& mat_p,
    const double alpha1, const double alpha2, bool is_debug) {
  if (is_debug) {
    logger().debug("calling update_tangent_points_on_tan_pls for mat_p {}...",
                   mat_p.tag);
    logger().debug("ref_fs_pairs_not_cross: {}", ref_fs_pairs_not_cross);
  }

  const int k = 2;  // k_ring sf_mesh face neighbors
  std::set<int> k_ring_fids;
  std::map<int, std::array<Vector3, 4>> faces_params;

  // store fids that are already stored
  // if another tangent plane found the same fids to update
  // then we delete duplicated tangent plane
  std::set<int> tangent_fids_now;
  std::vector<TangentPlane> tan_planes_new;

  for (auto& tan_pl : mat_p.tan_planes) {
    if (is_debug) {
      logger().debug("[Update TAN_POINT] before update:");
      tan_pl.print_info();
    }
    get_k_ring_neighbors_no_cross(sf_mesh, ref_fs_pairs_not_cross, tan_pl.fid,
                                  k, k_ring_fids, false /*is_debug*/);
    get_faces_params(sf_mesh, k_ring_fids, faces_params);
    // find tangent point that minimize the energy
    for (const auto& f_param_pair : faces_params) {
      int fid_tmp = f_param_pair.first;
      auto& f_param = f_param_pair.second;
      Vector3 pN_tmp;
      // Note that here we use alpha1=1, alpha2=0 to update tangent points
      double e_tmp = get_closest_tangent_point_opt(
          alpha1, alpha2, f_param[0], f_param[1], f_param[2], f_param[3],
          mat_p.pos, mat_p.sq_radius, pN_tmp, false /*is_debug*/);
      if (e_tmp < tan_pl.energy) {
        // update tangent point
        tan_pl.energy = e_tmp;
        tan_pl.fid = fid_tmp;
        tan_pl.points.clear();
        tan_pl.points.push_back(pN_tmp);
        tan_pl.normal = f_param[3];
        if (is_debug) {
          logger().debug(
              "[Update TAN_POINT] mat_p {} found lower energy {} with fid {}",
              mat_p.tag, e_tmp, fid_tmp);
        }
      }
    }  // for faces_params

    if (is_debug) {
      logger().debug("[Update TAN_POINT] after update:");
    }

    ////////////////
    // after update
    // 1. remove duplicated tangent planes
    if (tangent_fids_now.find(tan_pl.fid) != tangent_fids_now.end()) {
      if (is_debug) {
        logger().debug(
            "[Update TAN_POINT] mat_p {} with fid {} is stored, delete tan_pl",
            mat_p.tag, tan_pl.fid);
      }
      tan_pl.is_deleted = true;
      continue;
    }
    // 2. remove similar tangent planes
    for (const auto& tan_pl_new : tan_planes_new) {
      if (tan_pl_new.is_same_tan_pl(tan_pl)) {
        if (is_debug) {
          logger().debug(
              "[Update TAN_POINT] mat_p {} whose tangent plane already stored "
              "similar planes, skip storing",
              mat_p.tag);
        }
        tan_pl.is_deleted = true;
        break;
      }
    }
    if (!tan_pl.is_deleted) {
      tan_planes_new.push_back(tan_pl);
      tangent_fids_now.insert(tan_pl.fid);
    }
    if (is_debug) {
      tan_pl.print_info();
    }
  }  // for mat_p.tan_planes

  mat_p.tan_planes = tan_planes_new;
}

// return:
// true - break the loop
// false - keep looping
// is_good will tell us if it's good
//
// load tangent planes (pN, nN)
// load concave lines (c_pN, c_nN)
bool prepare_next_iterate(
    const int iteration_limit, const int num_itr, const double alpha1,
    const double alpha2, const double alpha3, const double break_threshold,
    const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper, MVertex& mat_p,
    double& sum_energy_over_sq_radius, std::vector<Vector3>& pN,
    std::vector<Vector3>& nN, std::vector<Vector3>& c_pN,
    std::vector<Vector3>& c_nN, bool& is_good, bool is_debug) {
  pN.clear();
  nN.clear();
  c_pN.clear();
  c_nN.clear();
  //////////////////////////
  // Break conditions: 2 parts
  if (num_itr != 0) {
    // check breaking condition
    // if num_itr != 0
    // fetch the nearest point q and qnormal
    Vector3 q;
    double _;
    int q_fid = aabb_wrapper.get_nearest_point_on_sf(mat_p.pos, q, _);
    const Vector3 qnormal = get_mesh_facet_normal(sf_mesh, q_fid);
    if (is_debug) {
      logger().debug(
          "[Iterate] mat_p {} num_itr {} has q_fid: {} with q: "
          "({},{},{}), qnormal: ({},{},{})",
          mat_p.tag, num_itr, q_fid, q[0], q[1], q[2], qnormal[0], qnormal[1],
          qnormal[2]);
    }
    // if q is too close to concave line
    // and qnormal is covered by mat_p.tan_cc_lines
    bool is_skip_add_new_tan_pl = false;
    double q_sq_dist_to_cc = aabb_wrapper.get_sq_dist_to_cc(q);
    if (q_sq_dist_to_cc < SQ_DIST_TO_CC &&
        mat_p.is_normal_covered_by_cc_lines(qnormal)) {
      is_skip_add_new_tan_pl = true;
      if (is_debug) {
        logger().debug(
            "[Iterate] mat_p {} num_itr: {} skip adding tan_pl of q_fid: {}, "
            "b/c "
            "too close to concave lines, q_sq_dist_to_cc: {}",
            mat_p.tag, num_itr, q_fid, q_sq_dist_to_cc);
      }
    }
    // create new tangent plane
    // and update existing tangent planes
    if (!is_skip_add_new_tan_pl) {
      TangentPlane tan_pl_q(qnormal);
      tan_pl_q.points.push_back(q);
      tan_pl_q.fid = q_fid;
      std::vector<TangentPlane> new_tan_planes;
      new_tan_planes.push_back(tan_pl_q);
      for (auto& tan_pl : mat_p.tan_planes) {
        if (tan_pl.is_same_normal(qnormal)) {
          tan_pl.is_deleted = true;
          if (is_debug) {
            logger().debug(
                "[Iterate] mat_p {} num_itr: {} remove tangent plane:",
                mat_p.tag, num_itr);
            tan_pl.print_info();
          }
        } else {
          new_tan_planes.push_back(tan_pl);
        }
      }
      if (is_debug)
        logger().debug(
            "[Iterate] mat_p {} num_itr {} added new tangent plane of q_fid "
            "{} "
            "q "
            "({},{},{}), qnormal ({},{},{})",
            mat_p.tag, num_itr, q_fid, q[0], q[1], q[2], qnormal[0], qnormal[1],
            qnormal[2]);
      mat_p.tan_planes = new_tan_planes;
    }  // if (!is_skip_add_new_tan_pl)

    //////////////////////////
    // Break condition part 1
    // iterate ends when the energy reaches minima
    double prev_energy_over_sq_radius = sum_energy_over_sq_radius;
    double sum_energy = mat_p.update_all_energy_values(alpha1, alpha2, alpha3);
    sum_energy_over_sq_radius = sum_energy / mat_p.sq_radius;

    if (is_debug) {
      logger().debug(
          "[Iterate] mat_p {} num_itr: {}, has sum_energy {}, "
          "sum_energy_over_sq_radius; {}",
          mat_p.tag, num_itr, sum_energy, sum_energy_over_sq_radius);
    }
    // this minima is small, good
    if (sum_energy_over_sq_radius < break_threshold) {
      if (is_debug)
        logger().debug(
            "[Iterate] sum_energy_over_sq_radius {} < break_threshold {}, "
            "minima is really small, good",
            sum_energy_over_sq_radius, break_threshold);

      is_good = true;
      return true;  // break
    }
    // change is too small, not good
    if (prev_energy_over_sq_radius != -1 &&
        std::abs(sum_energy_over_sq_radius - prev_energy_over_sq_radius) <=
            SCALAR_ZERO_N6) {
      if (is_debug)
        logger().debug(
            "[Iterate] sum_energy_over_sq_radius {} == "
            "prev_energy_over_sq_radius {}, break_threshold: {},"
            "change too small, not good",
            sum_energy_over_sq_radius, prev_energy_over_sq_radius,
            break_threshold);
      is_good = false;
      return true;  // break
    }

    //////////////////////////
    // Break condition part 2
    if (num_itr > iteration_limit || mat_p.sq_radius < SCALAR_ZERO_N6) {
      if (is_debug)
        logger().debug(
            "[Iterate] break, mat_p {} reaches iteration limits {} with sum "
            "energy: {}, mat_p.sq_radius {}",
            mat_p.tag, num_itr, sum_energy, mat_p.sq_radius);
      is_good = false;
      return true;  // break
    }
  }  // if (num_itr != 0)

  // For all iterations
  // prepare for next iteration
  for (const auto& tan_pl : mat_p.tan_planes) {
    pN.push_back(tan_pl.points[0]);
    nN.push_back(tan_pl.normal);
    if (is_debug) {
      logger().debug(
          "[Iterate] mat_p {} num_itr: {}, has tan_pl fid {} energy {}, "
          "reiterate",
          mat_p.tag, num_itr, tan_pl.fid, tan_pl.energy);
    }
  }  // for mat_p.tan_planes
  for (const auto& one_cc_line : mat_p.tan_cc_lines) {
    c_pN.push_back(one_cc_line.tan_point);
    c_nN.push_back(one_cc_line.normal);
    if (is_debug) {
      logger().debug(
          "[Iterate] mat_p {} num_itr: {}, has cc_line with energy {}, "
          "reiterate",
          mat_p.tag, num_itr, one_cc_line.energy);
    }
  }  // mat_p.tan_cc_lines
  is_good = false;
  return false;  // not break
}

// Given sphere (theta, r), find the closest tangent point pN
// and normal nN on concave segement (m1, m2) with adjacent normals (n1, n2)
// by solving an equation
double get_closest_concave_point(const double alpha3, const Vector3& m1,
                                 const Vector3& m2, const Vector3& n1,
                                 const Vector3& n2, const Vector3& theta,
                                 const double& sq_radius, Vector3& pN,
                                 Vector3& nN, bool is_debug) {
  double radius = std::sqrt(sq_radius);
  Vector3 m12 = m2 - m1;
  Vector3 vM = m12.normalized();
  double t = (theta.dot(vM) - m1.dot(vM)) / m12.dot(vM);
  double E = DBL_MAX;
  if (t >= 0 && t <= 1) {
    pN = m1 + t * m12;
    nN = (pN - theta).normalized();
    // nN must in range of (n1, n2)
    if (is_vector_within_range_of_two_vectors(n1, n2, nN)) {
      E = TangentConcaveLine::get_energy_value(theta, radius, alpha3, pN, nN);
      if (is_debug)
        logger().debug(
            "[CONCAVE_PNT] return pN ({},{},{}) with energy E {}, t: {}, nN: "
            "({},{},{})",
            pN[0], pN[1], pN[2], E, t, nN[0], nN[1], nN[2]);
      return E;
    }
  }
  // check boundary points/normals on segment
  // find the best combination with smallest energy
  // (pN, nN) pairs
  E = DBL_MAX;
  int pair_idx = -1;
  std::vector<std::array<Vector3, 2>> pN_nN_pairs{
      {{m1, nN}}, {{m2, nN}}, {{pN, n1}}, {{pN, n2}},
      {{m1, n1}}, {{m1, n2}}, {{m2, n1}}, {{m2, n2}}};
  for (int i = 0; i < pN_nN_pairs.size(); i++) {
    const auto& pair = pN_nN_pairs[i];
    double E_tmp = TangentConcaveLine::get_energy_value(theta, radius, alpha3,
                                                        pair[0], pair[1]);
    if (E_tmp < E) {
      pair_idx = i;
      E = E_tmp;
      pN = pair[0];
      nN = pair[1];
      // nN = (pN - theta).normalized();
    }
  }
  if (is_debug)
    logger().debug(
        "[CONCAVE_PNT] return pair_idx: {}, pN ({},{},{}) nN: ({},{},{}) with "
        "energy E {}, "
        "t: {}",
        pair_idx, pN[0], pN[1], pN[2], nN[0], nN[1], nN[2], E, t);

  return E;
}

void update_tangent_points_on_cc_lines(MVertex& mat_p, const double alpha3,
                                       bool is_debug) {
  if (is_debug && mat_p.tan_cc_lines.size() > 0)
    logger().debug("calling update_tangent_points_on_cc_lines for mat_p {}...",
                   mat_p.tag);

  for (auto& one_cc_line : mat_p.tan_cc_lines) {
    get_closest_concave_point(
        alpha3, one_cc_line.ref_vs_pos[0], one_cc_line.ref_vs_pos[1],
        one_cc_line.adj_ref_normals[0], one_cc_line.adj_ref_normals[1],
        mat_p.pos, mat_p.sq_radius, one_cc_line.tan_point, one_cc_line.normal,
        is_debug);
  }
}

/////////////////////////////////////////////////////////////////////////////////////
// Main Functions
/////////////////////////////////////////////////////////////////////////////////////
//
// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line
bool iterate_sphere(const GEO::Mesh& sf_mesh,
                    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
                    const AABBWrapper& aabb_wrapper, MVertex& mat_p,
                    bool is_debug, double alpha1, double alpha2, double alpha3,
                    const double break_threshold, const int iteration_limit) {
  if (is_debug) {
    logger().debug(
        "[Iterate] calling iterate_sphere for mat_p: {}, is_debug: {}",
        mat_p.tag, is_debug);
  }

  // update tangent planes first
  mat_p.update_tan_planes_by_sf_mesh(sf_mesh, aabb_wrapper);
  int num_itr = 0;
  // mat_p.store_sq_energy(aabb_wrapper);  // init
  double sum_energy_over_sq_radius = -1;
  std::vector<Vector3> pN, nN;      // for tangent planes
  std::vector<Vector3> c_pN, c_nN;  // for concave lines
  // alpha1 = 0.01;  // energy of distance to tangent point
  // alpha2 = 1;     // energy of distance to tangent plane
  // alpha3 = 1;     // energy of distance to concave line

  if (is_debug) {
    logger().debug(
        "[Iterate] before iteration, mat_p: {} has tangent planes and "
        "tan_cc_lines:",
        mat_p.tag);
    mat_p.print_tan_planes();
    mat_p.print_cc_lines();
  }

  bool is_good = false;
  while (true) {
    // update is_good
    // also update (pN, nN) and (c_pN, c_nN)
    bool is_break = prepare_next_iterate(
        iteration_limit, num_itr, alpha1, alpha2, alpha3, break_threshold,
        sf_mesh, aabb_wrapper, mat_p, sum_energy_over_sq_radius, pN, nN, c_pN,
        c_nN, is_good, is_debug);
    if (is_break) break;

    // Step 1: update sphere given tangent points
    Vector3 new_pos = Vector3::Zero();
    double new_sq_radius = 0.;
    bool is_success =
        update_sphere_given_plN_opt(alpha1, alpha2, pN, nN, alpha3, c_pN, c_nN,
                                    new_pos, new_sq_radius, is_debug);
    if (is_debug)
      logger().debug("update_sphere_given_plN_opt is_success: {}", is_success);
    if (!is_success || new_sq_radius < SCALAR_ZERO_N6) {
      //////////////////////////
      // Break condition part 4
      is_good = false;
      break;
    }
    mat_p.save_old_pos_radius();
    mat_p.pos = new_pos;
    mat_p.sq_radius = new_sq_radius;

    // Step 2: update tangent points given new sphere
    // 1. update tangent planes
    update_tangent_points_on_tan_pls(sf_mesh, ref_fs_pairs_not_cross, mat_p,
                                     alpha1, alpha2, is_debug /*is_debug*/);
    if (is_debug) {
      logger().debug(
          "[Update TAN_POINT] mat_p {} has tangent planes after update: ",
          mat_p.tag);
      mat_p.print_tan_planes();
    }
    // 2. update concave lines
    update_tangent_points_on_cc_lines(mat_p, alpha3, is_debug);
    num_itr++;
    // mat_p.store_sq_energy(aabb_wrapper);
  }  // while true

  if (is_debug)
    logger().debug("[Iterate] mat_p: {} iterated {} time, is_good: {}",
                   mat_p.tag, num_itr, is_good);
  if (is_good == false || num_itr == 0) {
    return false;
  }
  return true;
}

}  // namespace matfp