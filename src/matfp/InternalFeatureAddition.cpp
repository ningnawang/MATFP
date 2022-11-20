// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/InternalFeatureAddition.h"

#include "matfp/Common.h"
#include "matfp/InscribedSpheres.h"
#include "matfp/IterateSpheres.h"
#include "matfp/ShrinkSpheres.h"

namespace matfp {
/////////////////////////////////////////////////////////////////////////////////////
// Helper Functions
/////////////////////////////////////////////////////////////////////////////////////
void get_mat_A_and_B(MVertex& mat_p1, MVertex& mat_p2, MVertex*& mat_A,
                     MVertex*& mat_B) {
  // remove cc_part that cross concave lines
  // the connectivity around concave line is not reliable
  // we will check concave line later
  int p1_cc_parts = mat_p1.get_num_cc_parts_exclude_cc_lines();
  int p2_cc_parts = mat_p2.get_num_cc_parts_exclude_cc_lines();
  int p1_cc_line = mat_p1.tan_cc_lines.size();
  int p2_cc_line = mat_p2.tan_cc_lines.size();

  // mat_A always contains a cc_line if any
  if (p1_cc_line > 0 && p2_cc_line == 0) {
    mat_A = &mat_p1;
    mat_B = &mat_p2;
    return;
  }
  if (p2_cc_line > 0 && p1_cc_line == 0) {
    mat_A = &mat_p2;
    mat_B = &mat_p1;
    return;
  }

  // p1_cc_line == 0 && p2_cc_line == 0
  // or
  // p1_cc_line > 0 && p2_cc_line > 0
  // order here is to make sure duplicated edge
  // will not appear in visited_edges
  // now check cc_parts
  if (p1_cc_parts == p2_cc_parts) {  // ==
    if (mat_p1.valid_idx < mat_p2.valid_idx) {
      mat_A = &mat_p1;
      mat_B = &mat_p2;
    } else {
      mat_A = &mat_p2;
      mat_B = &mat_p1;
    }
    return;
  }

  if (p1_cc_parts < p2_cc_parts) {  // <
    mat_A = &mat_p1;
    mat_B = &mat_p2;
    return;
  }
  // >
  mat_A = &mat_p2;
  mat_B = &mat_p1;
  return;
}

// both A and B are common spheres T_N
// 1. when both T_2, then all two CCs must be adjacent,
//    if not, then create new T_N in add_a_new_sphere_v3()
// 2. when one T_2, other T_N, then at least 2 CCs must be adjacent,
//    if not, then create new T_2 in add_a_new_sphere_v3()
bool is_valid_if_any_T_N_sphere(bool is_faster, MVertex& mat_A, MVertex& mat_B,
                                const int adj_ccs, const int adj_tan_normals,
                                bool is_debug) {
  if (mat_A.type == SphereType::T_2 || mat_B.type == SphereType::T_2) {
    if ((is_faster && adj_tan_normals < 2) ||
        (!is_faster && (adj_ccs < 2 && adj_tan_normals < 2))) {
      if (is_debug) {
        logger().debug(
            "[INVALID CHECK] found invalid edge T_N ({},{}) has adj_ccs {}, "
            "adj_tan_normals {}",
            mat_A.tag, mat_B.tag, adj_ccs, adj_tan_normals);
      }
      return false;
    }
    return true;
  }
  return true;
}

bool is_valid_if_any_feature_sphere(bool is_faster, MVertex& mat_A,
                                    MVertex& mat_B, const int adj_ccs,
                                    const int adj_tan_normals, bool is_debug) {
  // if both are feature spheres then they should be adjacent
  if (mat_A.is_a_feature_sphere() && mat_B.is_a_feature_sphere()) {
    if (mat_A.se_adj_se.find(mat_B.tag) == mat_A.se_adj_se.end()) {
      if (is_debug)
        logger().debug(
            "[INVALID CHECK] found invalid edge feature ({},{}), not "
            "neighboring",
            mat_A.tag, mat_B.tag);
      return false;
    }
    return true;
  }

  if ((is_faster && adj_tan_normals < 2) ||
      (!is_faster && (adj_ccs < 2 || adj_tan_normals < 2))) {
    if (is_debug)
      logger().debug(
          "[INVALID CHECK] found invalid edge feature ({},{}), adj_ccs: "
          "{}, adj_tan_normals: {}",
          mat_A.tag, mat_B.tag, adj_ccs, adj_tan_normals);
    return false;
  }
  return true;
}

// For spheres contains cancave lines
//
// Each concave line has two adjacent normals
bool is_valid_if_any_T_2_c_sphere(bool is_faster, const MVertex& mat_A,
                                  const MVertex& mat_B, const int adj_ccs,
                                  const int adj_tan_normals, bool is_debug) {
  // if any T_2_c (= T_2_c) sphere
  // 1. if both T_2_c spheres, have 1 adjacent CC and 3 similar normals
  if (mat_A.type == SphereType::T_2_c && mat_B.type == SphereType::T_2_c) {
    // when adj_ccs = 1, the concave line could be completely different
    if ((is_faster && adj_tan_normals < 3) ||
        (!is_faster && (adj_ccs <= 1 && adj_tan_normals < 3))) {
      if (is_debug) {
        logger().debug(
            "[INVALID CHECK] found invalid concave edge ({},{}), adj_ccs: "
            "{}, adj_tan_normals: {}",
            mat_A.tag, mat_B.tag, adj_ccs, adj_tan_normals);
      }
      return false;
    }
    return true;
  }

  // 2. if one T_2_c sphere, then must adjacent to one CC
  if (mat_A.type == SphereType::T_2_c || mat_B.type == SphereType::T_2_c) {
    // if (adj_ccs < 1 || adj_tan_normals < 2) {
    if ((is_faster && adj_tan_normals < 2) ||
        (!is_faster && (adj_ccs < 1 || adj_tan_normals < 2))) {
      if (is_debug) {
        logger().debug(
            "[INVALID CHECK] found invalid concave edge ({},{}), adj_ccs: "
            "{}, adj_tan_normals: {}",
            mat_A.tag, mat_B.tag, adj_ccs, adj_tan_normals);
      }
      return false;
    }
  }
  return true;
}

///////////////////////////////////////////////////////
// Shared functions
///////////////////////////////////////////////////////

// order not matter, return max
int get_num_adj_cc_parts(const MVertex& mat_A, const MVertex& mat_B) {
  // Get number of adjacent CCs from A to B and B to A
  auto get_num_adj_cc_parts_A_to_B = [&](const MVertex& mat_A,
                                         const MVertex& mat_B) {
    int num_adj_cc_parts = 0;
    for (const auto& p1_cc : mat_A.cc_parts) {
      if (p1_cc.is_adj_to(mat_B.tag)) {
        num_adj_cc_parts++;
      }
    }
    return num_adj_cc_parts;
  };

  int num_ccs = get_num_adj_cc_parts_A_to_B(mat_A, mat_B);
  num_ccs = std::max(num_ccs, get_num_adj_cc_parts_A_to_B(mat_B, mat_A));
  return num_ccs;
}

// order not matter, return max
int get_num_common_tan_normals(const MVertex& mat_A, const MVertex& mat_B) {
  // Get number of similar tangent normals from A to B and B to A
  auto get_num_common_tan_normals_A_to_B =
      [&](std::vector<Vector3>& normals_A, std::vector<Vector3>& normals_B) {
        int num_common_tan_pls = 0;
        // yes, n^2
        for (const auto& A_normal : normals_A) {
          for (const auto& B_normal : normals_B) {
            if (is_vector_same_direction(A_normal, B_normal, esp_degree_20)) {
              num_common_tan_pls++;
              break;  // break normals_B
            }
          }
        }  // for normals_A
        return num_common_tan_pls;
      };

  std::vector<Vector3> normals_A, normals_B;
  mat_A.get_sphere_all_tangent_normals_includes_cc_lines(normals_A);
  mat_B.get_sphere_all_tangent_normals_includes_cc_lines(normals_B);

  int num_tan_normals = get_num_common_tan_normals_A_to_B(normals_A, normals_B);
  num_tan_normals = std::max(
      num_tan_normals, get_num_common_tan_normals_A_to_B(normals_B, normals_A));
  return num_tan_normals;
}

// the tangent planes of A and B should all coming from Shrink Spheres
// we aggregate them by checking the normal within certain threshold
void aggregate_tangent_planes(const MVertex& mat_A, const MVertex& mat_B,
                              MVertex& new_mat_p) {
  // copy A
  std::copy(mat_A.tan_planes.begin(), mat_A.tan_planes.end(),
            std::back_inserter(new_mat_p.tan_planes));
  // check B
  for (const auto& tan_pl_B : mat_B.tan_planes) {
    if (new_mat_p.is_normal_covered_then_merge(tan_pl_B)) {
      continue;
    } else {
      new_mat_p.tan_planes.push_back(tan_pl_B);
    }
  }
}

// TODO: dunno how to merge for now, save them all
void aggregate_concave_lines(const MVertex& mat_A, const MVertex& mat_B,
                             MVertex& new_mat_p) {
  // copy A
  std::copy(mat_A.tan_cc_lines.begin(), mat_A.tan_cc_lines.end(),
            std::back_inserter(new_mat_p.tan_cc_lines));
  // check B
  for (const auto& one_cc_line_B : mat_B.tan_cc_lines) {
    // if (new_mat_p.is_normal_covered_then_merge(tan_pl_B)) {
    //   continue;
    // } else {
    new_mat_p.tan_cc_lines.push_back(one_cc_line_B);
    // }
  }
}

// Find the common tangent planes, this is for one T_2 and another T_N
// consider adjacent tangent info of concave lines as well
void get_one_common_tangent_plane(const MVertex& mat_A, const MVertex& mat_B,
                                  MVertex& new_mat_p) {
  std::vector<std::array<Vector3, 2>> matA_tan_pairs, matB_tan_pairs;
  mat_A.get_sphere_all_tangent_pairs_includes_cc_lines(matA_tan_pairs);
  mat_B.get_sphere_all_tangent_pairs_includes_cc_lines(matB_tan_pairs);
  for (const auto& tan_pair_B : matB_tan_pairs) {
    const Vector3& B_point = tan_pair_B[0];
    const Vector3& B_normal = tan_pair_B[1];
    for (const auto& tan_pair_A : matA_tan_pairs) {
      const Vector3& A_point = tan_pair_A[0];
      const Vector3& A_normal = tan_pair_A[1];
      if (TangentPlane::is_same_normal(A_normal, B_normal)) {
        TangentPlane new_tan_pl(A_normal);
        new_tan_pl.points.push_back(A_point);
        new_tan_pl.points.push_back(B_point);
        new_mat_p.tan_planes.push_back(new_tan_pl);
        break;
      }
    }  // for matA_tan_pairs
  }    // for matB_tan_pairs
}

// create T_N sphere by calling iterate_sphere()
bool create_T_N_sphere(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper, const MVertex& mat_A, const MVertex& mat_B,
    MVertex& new_mat_p, bool is_debug) {
  if (is_debug)
    logger().debug(
        "[CREATE_TX] try adding T_N sphere for invalid edge A {} B {} ...",
        mat_A.tag, mat_B.tag);

  aggregate_tangent_planes(mat_A, mat_B, new_mat_p);
  new_mat_p.update_tan_planes_points();
  if (new_mat_p.tan_planes.size() < 2) {
    return false;
  }
  if (new_mat_p.tan_planes.size() == 2) {
    logger().error(
        "[CREATE_TX] aggregate tangent planes {}, not on internal feature, "
        "skip",
        new_mat_p.tan_planes.size());
    return true;
  }
  if (is_debug) new_mat_p.print_tan_planes();
  bool is_success =
      iterate_sphere(sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, new_mat_p,
                     is_debug /*is_debug*/);
  new_mat_p.type = SphereType::T_N;
  return is_success;
}

bool create_T_N_c_sphere(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper, const MVertex& mat_A, const MVertex& mat_B,
    MVertex& new_mat_p, bool is_debug) {
  if (is_debug)
    logger().debug(
        "[CREATE_TXc] try adding T_N_c sphere for invalid edge A {} B {} ...",
        mat_A.tag, mat_B.tag);
  aggregate_tangent_planes(mat_A, mat_B, new_mat_p);
  new_mat_p.update_tan_planes_points();
  aggregate_concave_lines(mat_A, mat_B, new_mat_p);
  if (is_debug) {
    new_mat_p.print_tan_planes();
    new_mat_p.print_cc_lines();
  }

  bool is_success =
      iterate_sphere(sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, new_mat_p,
                     is_debug /*is_debug*/);
  new_mat_p.type = SphereType::T_N_c;
  return is_success;
}

// create T_2 sphere by calling shrink_sphere_wrapper()
bool create_T_2_sphere(const GEO::Mesh& sf_mesh,
                       const AABBWrapper& aabb_wrapper, const MVertex& mat_A,
                       const MVertex& mat_B, MVertex& new_mat_p,
                       bool is_debug) {
  if (is_debug)
    logger().debug(
        "[CREATE_T2] try adding T_2 sphere for invalid edge A {} B {} ...",
        mat_A.tag, mat_B.tag);
  // create T_2 sphere
  get_one_common_tangent_plane(mat_A, mat_B, new_mat_p);
  new_mat_p.sq_radius = std::pow(INIT_RADIUS, 2);
  new_mat_p.update_tan_planes_by_sf_mesh(sf_mesh, aabb_wrapper);

  // will update to SphereType::T_2
  return shrink_sphere_wrapper(sf_mesh, aabb_wrapper, new_mat_p,
                               true /*is_setup_ss_param*/, true /*is_check_cc*/,
                               is_debug);
}

bool add_a_new_sphere_v3(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper, const MVertex& mat_A, const MVertex& mat_B,
    std::vector<MVertex>& all_medial_spheres,
    std::vector<MVertex>& sorted_partial_medial_spheres, bool is_shrink_only,
    bool is_debug) {
  if (is_debug)
    logger().debug(
        "[InternalAdd] adding one more for invalid edge ({}, {}) ...",
        mat_A.tag, mat_B.tag);

  // new mat
  int new_tag = all_medial_spheres.size();
  MVertex new_mat_p(new_tag, 0., 0., 0., 0.);

  bool is_success = false;
  if (is_shrink_only) {
    is_success = create_T_2_sphere(sf_mesh, aabb_wrapper, mat_A, mat_B,
                                   new_mat_p, is_debug);
  } else if (mat_A.is_a_concave_sphere() || mat_B.is_a_concave_sphere()) {
    is_success =
        create_T_N_c_sphere(sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper,
                            mat_A, mat_B, new_mat_p, is_debug);
  } else {
    is_success =
        create_T_N_sphere(sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, mat_A,
                          mat_B, new_mat_p, is_debug);
    // do nothing, no need to add more T_2 spheres
    if (is_success && new_mat_p.tan_planes.size() == 2) {
      return true;
    }
  }

  if (is_success) {
    // remember to update dilated radius
    new_mat_p.dilate_sphere_radius();
    new_mat_p.is_addition_int = AdditionType::A_FP_INT;
    new_mat_p.is_sphere_updated = true;
    all_medial_spheres.push_back(new_mat_p);
    sorted_partial_medial_spheres.push_back(new_mat_p);
    // if (is_debug)
    //   logger().debug(
    //       "[InternalAdd Created] invalid edge ({}, {}), new mat {}, cc_parts
    //       " "size: {}, " "tan_pl: {}, tan_cc_lines: {}, new pos ({},{},{})
    //       sq_radius {}, " "sq_radius_dilated: {}", mat_A.tag, mat_B.tag,
    //       new_tag, new_mat_p.cc_parts.size(), new_mat_p.tan_planes.size(),
    //       new_mat_p.tan_cc_lines.size(), new_mat_p.pos[0], new_mat_p.pos[1],
    //       new_mat_p.pos[2], new_mat_p.sq_radius,
    //       new_mat_p.sq_radius_dilated);
  } else {
    if (is_debug) {
      logger().debug(
          "[InternalAdd Uncreated] invalid edge ({}, {}), new mat {}, "
          "valid_idx: "
          "{}, cc_parts size: "
          "{}",
          mat_A.tag, mat_B.tag, new_tag, new_mat_p.valid_idx,
          new_mat_p.cc_parts.size());
      new_mat_p.print_info();
    }
  }  // if (is_success)

  return is_success;
}

/////////////////////////////////////////////////////////////////////////////////////
// Main Functions
/////////////////////////////////////////////////////////////////////////////////////
void check_invalid_mat_edges(
    NonManifoldMesh& mat, const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres,  // not const on purpose
    std::set<std::array<int, 2>>& invalid_mat_edges, bool is_faster,
    bool is_debug) {
  // NOTE: each edge in invalid_mat_edges
  //        must be sorted based on #CC
  // 	eg.   #CC(one_edge[0]) <= #CC(one_edge[1])
  //
  invalid_mat_edges.clear();
  std::set<std::array<int, 2>> visited_edges;
  for (const auto& bep : mat.edges) {
    int current_seed = (*bep.second).vertices_.first;
    int neig_seed = (*bep.second).vertices_.second;

    // check if mat edge is valid
    MVertex& mat_p1 = all_medial_spheres[valid_medial_spheres[current_seed]];
    MVertex& mat_p2 = all_medial_spheres[valid_medial_spheres[neig_seed]];
    // mat_A alwasy has less or equal cc_parts than mat_B
    MVertex* mat_A = nullptr;
    MVertex* mat_B = nullptr;  // copy
    get_mat_A_and_B(mat_p1, mat_p2, mat_A, mat_B);

    if (mat_A == nullptr || mat_B == nullptr) {
      continue;
    }

    std::array<int, 2> one_edge = {{mat_A->valid_idx, mat_B->valid_idx}};
    if (visited_edges.find(one_edge) != visited_edges.end()) continue;
    visited_edges.insert(one_edge);
    if (mat_A->type == SphereType::T_UNK && mat_B->type == SphereType::T_UNK) {
      continue;
    }

    int adj_ccs = get_num_adj_cc_parts(*mat_A, *mat_B);
    int adj_tan_normals = get_num_common_tan_normals(*mat_A, *mat_B);

    // when is_valid == true, we just continue the loop
    bool is_valid = true;
    // check if A/B is on sharp edge or corner
    if (mat_A->is_a_feature_sphere() || mat_B->is_a_feature_sphere()) {
      is_valid = is_valid_if_any_feature_sphere(
          is_faster, *mat_A, *mat_B, adj_ccs, adj_tan_normals, is_debug);
    }

    // // check if A/B on concave lines
    else if (mat_A->type == SphereType::T_2_c ||
             mat_B->type == SphereType::T_2_c) {
      is_valid = is_valid_if_any_T_2_c_sphere(
          is_faster, *mat_A, *mat_B, adj_ccs, adj_tan_normals, is_debug);
      if (!is_valid) invalid_mat_edges.insert(one_edge);
      // make sure we do not check after this
      continue;
    }
    ////////////
    // check if A&B are common types such as T2/T3
    //
    // This function will update mat.mat_intf_edges and mat.mat_intf_adjs
    else {
      is_valid = is_valid_if_any_T_N_sphere(is_faster, *mat_A, *mat_B, adj_ccs,
                                            adj_tan_normals, is_debug);
    }

    if (!is_valid) invalid_mat_edges.insert(one_edge);
    // logger().debug("[INVALID CHECK] found invalid edge tags ({},{}),
    // one_edge: {}", 	mat_A->tag, mat_B->tag, one_edge
    // );

  }  // for mat.edges

  logger().info("[INVALID CHECK] invalid mat edges size: {}",
                invalid_mat_edges.size());
}

void insert_internal_feature_spheres(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper,
    const std::set<std::array<int, 2>>& invalid_mat_edges,
    const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_debug) {
  if (is_debug)
    logger().debug("[INVALID DEBUG] looping {} invalid edges",
                   invalid_mat_edges.size());
  // make sure all new inserted spheres are unique
  std::vector<MVertex> sorted_partial_medial_spheres;

  // start looping invalid mat edges
  for (const auto& one_invalid_edge : invalid_mat_edges) {
    int sA = one_invalid_edge[0];
    int sB = one_invalid_edge[1];
    MVertex& mat_A = all_medial_spheres[valid_medial_spheres[sA]];
    MVertex& mat_B = all_medial_spheres[valid_medial_spheres[sB]];

    // 1. try to add T_N/T_N_c sphere by calling is_shrink_only = false
    bool is_success = add_a_new_sphere_v3(
        sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, mat_A, mat_B,
        all_medial_spheres, sorted_partial_medial_spheres,
        false /*is_shrink_only*/, is_debug);

    // 2. if failed, then add a T_2 sphere by calling is_shrink_only = true
    if (!is_success) {
      if (is_debug)
        logger().debug(
            "[INVALID] new mat {} creates T_N failed, create T_2 instead",
            all_medial_spheres.size());
      is_success = add_a_new_sphere_v3(
          sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, mat_A, mat_B,
          all_medial_spheres, sorted_partial_medial_spheres,
          true /*is_shrink_only*/, is_debug);
    }

    // 3. success, but sphere is too rare, add another T_2 just in case
    const MVertex& last_mat_p = all_medial_spheres.back();
    if (is_success && last_mat_p.tan_planes.size() >= 3) {
      if (is_debug)
        logger().error(
            "[INVALID] new mat {} of T_N/T_N_c will create another T_2",
            last_mat_p.tag);
      is_success = add_a_new_sphere_v3(
          sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, mat_A, mat_B,
          all_medial_spheres, sorted_partial_medial_spheres,
          true /*is_shrink_only*/, is_debug);
    }

  }  // for invalid_mat_edges

  logger().debug("[INVALID] added {} new spheres",
                 sorted_partial_medial_spheres.size());

  // Remove duplicates
  // sort spheres and delete duplicated spheres
  // if is_override = true: check all medial spehres, not only newly added
  remove_duplicated_medial_spheres(sorted_partial_medial_spheres,
                                   all_medial_spheres, false /*is_override*/);
}
}  // namespace matfp