// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "InscribedSpheres.h"

#include <assert.h>

#include <Eigen/Dense>
#include <algorithm>
#include <random>

#include "matfp/Common.h"
#include "matfp/Logger.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////////
// Common functions used by UpdateSpheres.h and IterateSpheres.h
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////
// Helper functions
// Non-feature tets
// these function are called after DTs
TET_TYPE get_nonfeature_tet_type(
    const std::vector<int>& tet_vs_tags,
    const std::vector<Vector3>& tet_vs_normals,
    std::vector<std::vector<int>>&
        tet_vs_bucket, /*group vs by normals, store 1-4 only*/
    bool is_debug = false) {
  // std::vector<std::vector<int>> tet_vs_bucket;
  std::unordered_set<int> visited;
  for (int i = 0; i < 4; i++) {
    if (visited.find(i) != visited.end()) continue;
    std::vector<int> sub_bucket;
    sub_bucket.push_back(i);
    visited.insert(i);
    Vector3 n_i = tet_vs_normals[i];  // copy
    for (int j = i + 1; j < 4; j++) {
      if (visited.find(j) != visited.end()) continue;
      Vector3 n_j = tet_vs_normals[j];  // copy
      if (is_vector_same_direction(n_i, n_j, esp_degree_10)) {
        sub_bucket.push_back(j);
        visited.insert(j);
      }
    }
    tet_vs_bucket.push_back(sub_bucket);
  }

  std::vector<int> b_sizes;
  for (const auto& b : tet_vs_bucket) {
    b_sizes.push_back(b.size());
  }
  std::sort(b_sizes.begin(), b_sizes.end(),
            std::greater<int>());  // descending order
  assert(b_sizes.size() > 0);

  if (is_debug) {
    logger().error(
        "non-feature tet tet_vs_tags: {} tet_vs_bucket {}, b_sizes {}",
        tet_vs_tags, tet_vs_bucket, b_sizes);
    for (int i = 0; i < tet_vs_normals.size(); i++) {
      logger().error("normal {}: ({},{},{})", i, tet_vs_normals[i][0],
                     tet_vs_normals[i][1], tet_vs_normals[i][2]);
    }
  }

  if (b_sizes[0] == 4)
    return TET_TYPE::T4_0;
  else if (b_sizes[0] == 3)
    return TET_TYPE::T3_1;
  else if (b_sizes[0] == 2) {
    assert(b_sizes.size() > 1);
    // check next
    if (b_sizes[1] == 2)
      return TET_TYPE::T2_2;
    else if (b_sizes[1] == 1)
      return TET_TYPE::T2_1_1;
    else {
      logger().error("b_sizes: {}", b_sizes);
      log_and_throw("ERROR: do not know how to classify this tet");
    }
  } else if (b_sizes[0] == 1) {
    if (b_sizes.size() == 4) {
      return TET_TYPE::T1_1_1_1;
    } else {
      logger().error("b_sizes: {}", b_sizes);
      log_and_throw("ERROR: do not know how to classify this tet");
    }
  } else {
    logger().error("b_sizes: {}", b_sizes);
    log_and_throw("ERROR: do not know how to classify this tet");
  }

  return TET_TYPE::T4_0;
}

void sort_tet_vs_and_normals(const std::vector<std::vector<int>>& tet_vs_bucket,
                             const std::vector<Vector3>& tet_vs_pos,
                             const std::vector<Vector3>& tet_vs_normals,
                             std::vector<Vector3>& tet_pos_sorted,
                             std::vector<Vector3>& tet_normal_sorted) {
  tet_pos_sorted.clear();
  tet_normal_sorted.clear();

  // add 3 or 2 first
  for (const auto& sub_bucket : tet_vs_bucket) {
    if (sub_bucket.size() == 3 || sub_bucket.size() == 2) {
      for (const auto& i : sub_bucket) {
        tet_pos_sorted.push_back(tet_vs_pos[i]);
        tet_normal_sorted.push_back(tet_vs_normals[i]);
      }
    }
  }

  // add 1
  for (const auto& sub_bucket : tet_vs_bucket) {
    if (sub_bucket.size() == 1) {
      for (const auto& i : sub_bucket) {
        tet_pos_sorted.push_back(tet_vs_pos[i]);
        tet_normal_sorted.push_back(tet_vs_normals[i]);
      }
    }
  }

  if (tet_pos_sorted.size() != 4 || tet_normal_sorted.size() != 4) {
    logger().debug(
        "tet_pos_sorted size {} != 4, OR tet_normal_sorted.size {} != 4",
        tet_pos_sorted.size(), tet_normal_sorted.size());
    log_and_throw("ERROR");
  }
}

void collect_tangent_planes_for_DT(
    const TET_TYPE tet_type, const std::vector<Vector3>& tet_vs_pos,
    const std::vector<Vector3>& tet_vs_normals,
    const std::vector<std::vector<int>>& tet_vs_bucket,
    std::vector<TangentPlane>& tan_elements,
    std::vector<int>& selected_tan_pls) {
  // eg. 3:1, then sorted as ABC, D.
  // eg. 2:2, then sorted as AB, CD
  // eg. 2:1:1, then sorted as AB, C, D
  // eg. 1:1:1:1, then sorted as A, B, C, D
  std::vector<Vector3> tet_pos_sorted;     // pA, pB, pC, pD
  std::vector<Vector3> tet_normal_sorted;  // nA, nB, nC, nD

  tan_elements.clear();
  selected_tan_pls.clear();
  sort_tet_vs_and_normals(tet_vs_bucket, tet_vs_pos, tet_vs_normals,
                          tet_pos_sorted, tet_normal_sorted);

  if (tet_type == TET_TYPE::T3_1) {
    TangentPlane tan_e1(tet_normal_sorted[0]);
    Vector3 mid =
        1. / 3. * (tet_pos_sorted[0] + tet_pos_sorted[1] + tet_pos_sorted[2]);
    tan_e1.push_new_point(mid);
    TangentPlane tan_e2(tet_normal_sorted[3]);
    tan_e2.push_new_point(tet_pos_sorted[3]);
    tan_elements.push_back(tan_e1);
    tan_elements.push_back(tan_e2);
    selected_tan_pls.push_back(0);
    selected_tan_pls.push_back(1);
    return;
  } else if (tet_type == TET_TYPE::T2_2) {
    TangentPlane tan_e1(tet_normal_sorted[0]);
    Vector3 mid1 = 1. / 2. * (tet_pos_sorted[0] + tet_pos_sorted[1]);
    tan_e1.push_new_point(mid1);
    TangentPlane tan_e2(tet_normal_sorted[2]);
    Vector3 mid2 = 1. / 2. * (tet_pos_sorted[2] + tet_pos_sorted[3]);
    tan_e2.push_new_point(mid2);
    tan_elements.push_back(tan_e1);
    tan_elements.push_back(tan_e2);
    selected_tan_pls.push_back(0);
    selected_tan_pls.push_back(1);
  } else if (tet_type == TET_TYPE::T2_1_1) {
    TangentPlane tan_e1(tet_normal_sorted[0]);
    Vector3 mid = 1. / 2. * (tet_pos_sorted[0] + tet_pos_sorted[1]);
    tan_e1.push_new_point(mid);
    TangentPlane tan_e2(tet_normal_sorted[2]);
    tan_e2.push_new_point(tet_pos_sorted[2]);
    TangentPlane tan_e3(tet_normal_sorted[3]);
    tan_e3.push_new_point(tet_pos_sorted[3]);
    tan_elements.push_back(tan_e1);
    tan_elements.push_back(tan_e2);
    tan_elements.push_back(tan_e3);
    selected_tan_pls.push_back(0);
    selected_tan_pls.push_back(1);
    selected_tan_pls.push_back(2);
  } else if (tet_type == TET_TYPE::T1_1_1_1) {
    TangentPlane tan_e1(tet_normal_sorted[0]);
    tan_e1.push_new_point(tet_pos_sorted[0]);
    TangentPlane tan_e2(tet_normal_sorted[1]);
    tan_e2.push_new_point(tet_pos_sorted[1]);
    TangentPlane tan_e3(tet_normal_sorted[2]);
    tan_e3.push_new_point(tet_pos_sorted[2]);
    TangentPlane tan_e4(tet_normal_sorted[3]);
    tan_e4.push_new_point(tet_pos_sorted[3]);
    tan_elements.push_back(tan_e1);
    tan_elements.push_back(tan_e2);
    tan_elements.push_back(tan_e3);
    tan_elements.push_back(tan_e4);
    selected_tan_pls.push_back(0);
    selected_tan_pls.push_back(1);
    selected_tan_pls.push_back(2);
    selected_tan_pls.push_back(3);
  }
}

//////////////////////////////////////////////////////
// Main function:
// Update sphere right after DT
void collect_and_update_DT_feature_tets(
    const DelaunayTriangulation& dt, const std::vector<double>& seed_points,
    const std::vector<Vector3>& seed_normals_v2,
    std::vector<MVertex>& all_medial_spheres) {
  logger().debug("start updating tet spheres ...");

  // tet_vs_tags = v_feature_tags + v_non_feature_tags
  std::vector<int> tet_vs_tags(4, -1);  // store tags of all vertices
  std::vector<Vector3> tet_vs_pos(4);
  std::vector<Vector3> tet_vs_normals(4);
  std::vector<int> v_feature_tags;      // store tags of feature vertices
  std::vector<int> v_non_feature_tags;  // store tags of non-feature vertices

  ////////////
  // start looping all tetrahedrons
  ////////////
  for (Finite_cells_iterator_dt fci = dt.finite_cells_begin();
       fci != dt.finite_cells_end(); fci++) {
    if (fci->info().is_outside) continue;

    int cell_tag = fci->info().tag;
    MVertex& mat_p = all_medial_spheres[cell_tag];
    if (mat_p.is_outside)
      // log_and_throw("ERROR: mat_p cannot be outside while updating in DT");
      continue;

    int nb_feature_vs = 0;
    v_feature_tags.clear();
    v_non_feature_tags.clear();
    bool is_contain_corner = false;
    bool is_contain_c_mesh_vs = false;
    bool is_outside = false;
    for (int i = 0; i < 4; i++) {
      tet_vs_tags[i] = fci->vertex(i)->info().tag;
      tet_vs_pos[i] = fci->vertex(i)->info().pos;
      if (tet_vs_tags[i] >= seed_normals_v2.size()) {  // TODO: check why
        logger().error(
            "tet sphere {} contains vertex {} that is >= seed_normals_v2.size, "
            "dunno why"
            "{}",
            mat_p.tag, tet_vs_tags[i], seed_normals_v2.size());
        mat_p.is_outside = true;
        is_outside = true;
        break;
      }
      tet_vs_normals[i] = seed_normals_v2[tet_vs_tags[i]];

      if (fci->vertex(i)->info().is_on_feature) {
        v_feature_tags.push_back(tet_vs_tags[i]);
        nb_feature_vs++;
      } else {
        v_non_feature_tags.push_back(tet_vs_tags[i]);
      }
    }
    if (is_outside) {
      continue;
    }

    // if (!is_contain_corner)
    // 	continue;
    if (nb_feature_vs > 4)
      log_and_throw("ERROR: tetra cannot have > 4 vertices on feature");

    ////////////////////////////////////////////////
    // updating non-feature tets
    if (nb_feature_vs <= 0) {
      /*group tet vs by normals, store 1-4 only*/
      std::vector<std::vector<int>> tet_vs_bucket;
      bool is_debug = false;
      TET_TYPE tet_type = get_nonfeature_tet_type(tet_vs_tags, tet_vs_normals,
                                                  tet_vs_bucket, is_debug);
      if (is_debug) {
        logger().error("cell_tag {} tet_type: {}", cell_tag,
                       static_cast<int>(tet_type));
      }
      // 4:0 is sliver, we do not know how to handle
      if (tet_type == TET_TYPE::T4_0) {
        mat_p.is_deleted = true;
        continue;
      }
      collect_tangent_planes_for_DT(tet_type, tet_vs_pos, tet_vs_normals,
                                    tet_vs_bucket, mat_p.tan_planes,
                                    mat_p.selected_tan_pls);
      mat_p.is_to_be_updated = true;
      continue;
    }

    // Skip feature tets
    // deprecating rest of the code!!! just skip
    continue;
  }  // for each tetra(dt cell)

  logger().debug("done updating tet spheres ");
}

///////////////////////////////////////////////////
// Main functions
// used by UpdateSpheres.h and IterateSpheres.h
///////////////////////////////////////////////////
void update_intf_adjacency(const int cur_tag,
                           std::vector<MVertex>& all_medial_spheres) {
  // fetch from all_medial_spheres again
  // which has been updated
  // otherwise sorted_partial_medial_spheres is not updated
  const auto& mat_curr = all_medial_spheres[cur_tag];
  if (mat_curr.intf_adj_intf.size() != 2) {
    for (const int& inf_adj : mat_curr.intf_adj_intf) {
      auto& mat_adj = all_medial_spheres[inf_adj];
      mat_adj.intf_adj_intf.erase(mat_curr.tag);
    }
    return;
  }

  std::vector<int> intf_adj_intf_vec;
  intf_adj_intf_vec.insert(intf_adj_intf_vec.end(),
                           mat_curr.intf_adj_intf.begin(),
                           mat_curr.intf_adj_intf.end());
  int adj_tag_1 = intf_adj_intf_vec[0];
  int adj_tag_2 = intf_adj_intf_vec[1];
  auto& mat_adj_1 = all_medial_spheres[adj_tag_1];
  auto& mat_adj_2 = all_medial_spheres[adj_tag_2];
  mat_adj_1.intf_adj_intf.erase(mat_curr.tag);
  mat_adj_2.intf_adj_intf.erase(mat_curr.tag);

  mat_adj_1.intf_adj_intf.insert(adj_tag_2);
  mat_adj_2.intf_adj_intf.insert(adj_tag_1);

  logger().debug(
      "[INTF renew] mat_curr {} is deleted, update intf_adj_intf: {}",
      mat_curr.tag, mat_curr.intf_adj_intf);
  logger().debug("[INTF renew] mat_adj_1 {} has intf_adj_intf: {}",
                 mat_adj_1.tag, mat_adj_1.intf_adj_intf);
  logger().debug("[INTF renew] mat_adj_2 {} has intf_adj_intf: {}",
                 mat_adj_2.tag, mat_adj_2.intf_adj_intf);
}

// will sort sorted_partial_medial_spheres in place
int remove_duplicated_medial_spheres(
    std::vector<MVertex>& sorted_partial_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_override) {
  if (is_override || sorted_partial_medial_spheres.empty()) {
    sorted_partial_medial_spheres.clear();
    sorted_partial_medial_spheres = all_medial_spheres;
  }
  // make sure all pos of mat_p are accessable before sorting
  for (auto& mat_p : all_medial_spheres) {
    mat_p.ensure_valid();
  }

  // Remove duplicates
  // sort spheres and delete duplicated spheres
  logger().debug("[Duplicate] start sorting given medial spheres ...");
  std::sort(sorted_partial_medial_spheres.begin(),
            sorted_partial_medial_spheres.end());
  int curr_idx = 0;  // init
  std::set<int> deleted_set;
  for (int n_idx = 1; n_idx < sorted_partial_medial_spheres.size(); n_idx++) {
    MVertex& mat_curr = sorted_partial_medial_spheres[curr_idx];
    MVertex& mat_next = sorted_partial_medial_spheres[n_idx];
    // sanity check
    if (mat_curr != all_medial_spheres[mat_curr.tag]) {
      mat_curr.print_info();
      all_medial_spheres[mat_curr.tag].print_info();
      log_and_throw("[Duplicate] sorted_mat != all_medial_spheres");
    }
    // if any already deleted
    if (mat_curr.is_deleted) {
      curr_idx = n_idx;
      continue;
    }
    if (mat_next.is_deleted) {
      continue;
    }
    // mat_next would have higher type than mat_curr,
    // so we delete mat_curr
    if (mat_curr == mat_next) {
      if (mat_curr.intf_adj_intf.size() <= 2) {
        // this is not a internal junction, can be deleted
        // mark delete
        mat_curr.is_deleted = true;  // this is just a copy
        all_medial_spheres[mat_curr.tag].is_deleted =
            true;  // this is the real deletion
        all_medial_spheres[mat_curr.tag].is_deleted_int = DeletionType::D_DUP;
        deleted_set.insert(mat_curr.tag);
        // logger().debug("[Delete] delete mat_curr {}, mat_next {}",
        // mat_curr.tag, mat_next.tag);

        // update internal feature adjacency
        update_intf_adjacency(mat_curr.tag, all_medial_spheres);
      } else {
        logger().debug(
            "[Delete] mat_curr {} is an internal junction, cannot be deletd, "
            "intf_adj_intf: {}, mat_next {}",
            mat_curr.tag, mat_curr.intf_adj_intf, mat_next.tag);
      }
    }
    curr_idx = n_idx;
  }
  logger().debug("sorted_partial_medial_spheres: {}, all_medial_spheres: {}",
                 sorted_partial_medial_spheres.size(),
                 all_medial_spheres.size());
  logger().info("Removed {} duplicated MAT spheres", deleted_set.size());
  return deleted_set.size();
}

bool is_internal_medial_feature(const SphereType& type) {
  // not a feature
  if (type <= SphereType::T_2) {
    return false;
  }
  // skip if only contain 1 cc_line
  if (type == SphereType::T_2_c) {
    return false;
  }
  return true;
}

}  // namespace matfp
