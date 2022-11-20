// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/MedialSpheresProcessor.h"

#include <memory>
#include <queue>

#include "matfp/InscribedSpheres.h"
#include "matfp/InternalFeatureAddition.h"
#include "matfp/IterateSpheres.h"
#include "matfp/ThreeDimensionalShape.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////
// Helper Functions: Sharp Edges
////////////////////////////////////////////////////////////////////////////////
void reload_se_spheres_adjs(
    const std::set<std::array<int, 2>>& se_spheres,
    std::map<int, std::unordered_set<int>>& se_spheres_adjs) {
  se_spheres_adjs.clear();
  // save s_edges in another way
  for (auto const& se_pair : se_spheres) {
    se_spheres_adjs[se_pair[0]].insert(se_pair[1]);
    se_spheres_adjs[se_pair[1]].insert(se_pair[0]);
  }
}

// this function is matching se_spheres now
// remember to call reload_se_spheres_adjs ahead!!
void reload_se_kd_tree(
    const std::vector<MVertex>& all_medial_spheres,
    const std::set<std::array<int, 2>>& se_spheres,
    const std::map<int, std::unordered_set<int>>& se_spheres_adjs,
    GEO::NearestNeighborSearch_var& se_kd_tree,
    std::vector<double>& se_kd_points,
    std::map<int, int>& se_kd_tree_idx_to_se_tag, bool is_debug) {
  if (se_kd_tree == nullptr) {
    se_kd_tree = GEO::NearestNeighborSearch::create(3);
    logger().debug("created se_kd_tree is null? {}", se_kd_tree == nullptr);
  }
  // se_kd_tree.reset(); // reset will make it nullptr
  se_kd_points.clear();
  // idx from se_kd_points to se_spheres_adjs keys
  se_kd_tree_idx_to_se_tag.clear();
  for (const auto se_pair : se_spheres_adjs) {
    int se_tag = se_pair.first;
    int se_p_kd_idx = se_kd_points.size() / 3;
    const Vector3 pos = all_medial_spheres.at(se_tag).pos;
    for (int i = 0; i < 3; i++) {
      se_kd_points.push_back(pos[i]);
    }
    se_kd_tree_idx_to_se_tag[se_p_kd_idx] = se_tag;
  }
  se_kd_tree->set_points(se_kd_points.size() / 3, se_kd_points.data());
  if (is_debug)
    logger().debug("[se_kd_tree reload] tree has nb_points: {}",
                   se_kd_tree->nb_points());
}

void fetch_k_nearest_sharp_edges(
    const size_t K_nb_neigh,  // k nearest neighbors
    const std::vector<MVertex>& all_medial_spheres,
    const std::map<int, std::unordered_set<int>>& se_spheres_adjs,
    const GEO::NearestNeighborSearch_var& se_kd_tree,
    const std::map<int, int>& se_kd_tree_idx_to_se_tag, const Vector3& center,
    std::set<std::array<int, 2>>& k_nearest_se, bool is_debug) {
  if (is_debug) logger().debug("calling fetch_k_nearest_sharp_edges ...");
  k_nearest_se.clear();
  // if no sharp edge, then skip
  if (se_kd_tree->nb_points() < 1) return;
  // fetch k nearest points on sharp edge
  std::vector<GEO::index_t> k_neighbors(K_nb_neigh);
  std::vector<double> k_neighbors_sq_dist(K_nb_neigh);
  std::vector<double> pd = {center[0], center[1], center[2]};
  se_kd_tree->get_nearest_neighbors(K_nb_neigh, pd.data(), k_neighbors.data(),
                                    k_neighbors_sq_dist.data());
  // logger().debug("kd-tree found {} nearest se points: {}", K_nb_neigh,
  // k_neighbors);
  std::vector<bool> is_se_p_checked(se_spheres_adjs.size(), false);
  for (const GEO::index_t kn_se_id : k_neighbors) {
    const int se_tag = se_kd_tree_idx_to_se_tag.at(kn_se_id);
    if (is_debug)
      logger().debug("kd-tree found nearest se point {}, mapped to idx {}",
                     kn_se_id, se_tag);
    if (se_spheres_adjs.find(se_tag) == se_spheres_adjs.end()) {
      logger().error("Error: se_tag {} not in se_spheres_adjs, kn_se_id: {}",
                     se_tag, kn_se_id);
      log_and_throw("ERROR: We cannot find sharp edge in se_spheres_adjs");
    }
    const std::unordered_set<int>& adj_se = se_spheres_adjs.at(se_tag);
    for (const auto& adj_se_i : adj_se) {
      // each se_tag might have more than 2 adj_se_i
      // store one sharp edge
      std::array<int, 2> one_edge = {{se_tag, adj_se_i}};
      std::sort(one_edge.begin(), one_edge.end());
      k_nearest_se.insert(one_edge);
    }
  }  // for k_neighbors
}

// AB is a feature edge indicating by two points A and B
// X is a mat point with sq_radius
// it would return if this mat point should be delete or not
//
// dilation won't impact this, because we just add more if r_term is smaller
bool is_rpd_break_one_sharp_edge(const Vector3& A, const Vector3& B,
                                 const Vector3& X, const double& sq_radius,
                                 bool is_debug = false) {
  double l_term = -A.dot(B);
  double r_term = X.dot(X) - X.dot(A + B) - sq_radius;
  if (is_debug) logger().debug("l_term {}, r_term {}", l_term, r_term);

  // if l_term < r_term, then no need to delete
  // else, delete
  double robust = SCALAR_ZERO_3;
  if (l_term - r_term + robust < 0.0) return false;
  return true;
};

// TODO:
// can we only detect neighboring sharp edges, not all?
bool is_delete_mat_given_all_sharp_edges(
    const std::vector<MVertex>& all_medial_spheres,
    const std::set<std::array<int, 2>>& se_spheres, const Vector3& X,
    const double& sq_radius) {
  for (auto e : se_spheres) {
    Vector3 v0 = all_medial_spheres.at(e[0]).pos;
    Vector3 v1 = all_medial_spheres.at(e[1]).pos;
    // logger().debug("checking mat ({},{},{}), processing feature edge v0:
    // {}-({},{},{}), v1: {}-({},{},{})",
    //     X[0], X[1], X[2],
    //     e[0], v0[0], v0[1], v0[2], e[1], v1[0], v1[1], v1[2]
    // );
    if (is_rpd_break_one_sharp_edge(v0, v1, X, sq_radius)) return true;
  }
  return false;
}

// Remove medial spheres if too close to sharp edges
// this function would only remove, will not add
//
// TODO
// can we only detect neighboring sharp edges, not all?
void remove_medial_spheres_too_close_to_sharp_edges(
    const std::set<std::array<int, 2>>& se_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_using_dilated_radius,
    bool is_debug) {
  logger().debug("removing spheres that too close to sharp edges ...");

  ///////////
  // checking all mat points
  std::set<int> deleted_spheres;
  for (int i = 0; i < all_medial_spheres.size(); i++) {
    MVertex& mat_p = all_medial_spheres[i];
    if (mat_p.is_outside) continue;

    // lock feature mat pts
    if (mat_p.is_a_feature_sphere()) {
      mat_p.is_locking = true;
      continue;
    }
    // if the mat is already deleted, continue
    if (mat_p.is_deleted) {
      continue;
    }

    // filter non-feature INSIDE mat pts
    // TODO: use k_nearest_se instead
    double sq_radius = mat_p.sq_radius;
    if (is_using_dilated_radius) {
      mat_p.dilate_sphere_radius();
      sq_radius = mat_p.sq_radius_dilated;
    }
    if (is_delete_mat_given_all_sharp_edges(all_medial_spheres, se_spheres,
                                            mat_p.pos, sq_radius)) {
      mat_p.is_deleted_for_se = true;
      mat_p.is_deleted = true;
      mat_p.is_deleted_int = DeletionType::D_FP_EXT;
      deleted_spheres.insert(mat_p.tag);
      continue;
    }
  }  // for all_medial_spheres done
  logger().debug("deleted #{}/{} medial spheres", deleted_spheres.size(),
                 all_medial_spheres.size());
}

bool check_and_add_to_se_recursively(
    const std::array<MVertex, 2> one_se_spheres,  // in tag
    const MVertex& target_sphere, std::vector<MVertex>& new_feature_spheres,
    std::set<std::array<int, 2>>& se_spheres_new, int& latest_tag,
    bool is_using_dilated_radius, bool is_debug) {
  if (is_debug) {
    logger().debug("[se_rec] calling check_and_add_to_se_recursively ...");
    logger().debug("[se_rec] checking se: ({},{})", one_se_spheres[0].tag,
                   one_se_spheres[1].tag);
  }
  Vector3 v0 = one_se_spheres[0].pos;
  Vector3 v1 = one_se_spheres[1].pos;
  std::array<int, 2> se_copy = {{one_se_spheres[0].tag, one_se_spheres[1].tag}};
  std::sort(se_copy.begin(), se_copy.end());
  // if the sharp edge is degenerated to a point
  // this actually means no matter how many points added in between,
  // target_sphere would break the connectivity of sharp edge (v0, v1);
  // TODO: do something other than return
  if ((v0 - v1).norm() < SCALAR_ZERO_N3) {
    if (is_debug)
      logger().error(
          "[se_rec] target_sphere {} cannot add new point to se_copy {}",
          target_sphere.tag, se_copy);
    // save old sharp edges
    se_spheres_new.insert(se_copy);
    return false;
  }

  // aggregates normals for new feature sphere
  std::vector<Vector3> se_normals;
  std::vector<int> se_ref_fids;  // matching GEO::Mesh
  for (const auto& mat_p : one_se_spheres) {
    for (const auto& tan_pl : mat_p.tan_planes) {
      se_normals.push_back(tan_pl.normal);
      se_ref_fids.push_back(tan_pl.fid);
    }
  }
  double sq_radius = target_sphere.sq_radius;
  if (is_using_dilated_radius) {
    /* use dilated radius!!!!*/
    sq_radius = target_sphere.sq_radius_dilated;
  }
  bool is_break =
      is_rpd_break_one_sharp_edge(v0, v1, target_sphere.pos, sq_radius);
  if (is_break) {
    // density is not enough
    // add midpoint to sample
    Vector3 mid = 1. / 2. * (v0 + v1);
    MVertex mid_sphere(latest_tag++, mid[0], mid[1], mid[2],
                       SCALAR_FEATURE_SQ_RADIUS, SphereType::T_1_N);
    mid_sphere.is_on_s_edge = true;
    mid_sphere.type == SphereType::T_1_N;
    mid_sphere.dilate_sphere_radius();
    mid_sphere.add_tan_planes_for_feature_sphere(se_normals, se_ref_fids);
    new_feature_spheres.push_back(mid_sphere);
    if (is_debug)
      logger().debug("[se_rec] bad, add new mid_sphere: {}", mid_sphere.tag);
    // left part
    bool is_good = check_and_add_to_se_recursively(
        {{one_se_spheres[0], mid_sphere}}, target_sphere, new_feature_spheres,
        se_spheres_new, latest_tag, is_using_dilated_radius, is_debug);
    if (!is_good) return false;
    // right part
    is_good = check_and_add_to_se_recursively(
        {{mid_sphere, one_se_spheres[1]}}, target_sphere, new_feature_spheres,
        se_spheres_new, latest_tag, is_using_dilated_radius, is_debug);
    if (!is_good) return false;
  } else {
    // save/update sharp edges
    if (is_debug) logger().debug("[se_rec] good, save new se_copy {}", se_copy);
    se_spheres_new.insert(se_copy);
    return true;
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
// Helper Functions: Corners
////////////////////////////////////////////////////////////////////////////////

// find small area close to corner to perform operations
void load_corners_min_length(const std::vector<SpecialEdge>& sp_edges,
                             SpecialCorner& sp_corner,
                             std::vector<MVertex>& all_medial_spheres) {
  // loop both sharp edges and concave edges
  const auto& adj_edges_sorted = sp_corner.adj_edges_sorted;
  int num_adj_edges = adj_edges_sorted.size();
  for (int cur_idx = 0; cur_idx < num_adj_edges; cur_idx++) {
    auto& cur_spe = sp_edges.at(adj_edges_sorted.at(cur_idx));
    // min distance between corner and centroids
    for (int i = 0; i < 2; i++) {
      double len_tmp = get_distance_between_two_vectors(
          sp_corner.pos, cur_spe.adj_ref_tan_points[i]);
      sp_corner.min_corner_len = std::min(sp_corner.min_corner_len, len_tmp);
    }
  }
  sp_corner.min_corner_len /= 3;  // make this distance smaller

  auto& mat_corner = all_medial_spheres.at(sp_corner.mat_tag);
  mat_corner.min_corner_len = sp_corner.min_corner_len;
}

// Get all sheet nodes defined by special edges
// we loop SE first to define sheets with 2 tangent planes
// then loop CCE to add sheets with concave lines
// store a queue of SpecialEdge for tracing sheet nodes extends from seam
void load_sheet_nodes_sorted_for_sp_corner(
    const SpecialCorner& sp_corner,
    const std::vector<TangentConcaveLine> tan_cc_lines,
    std::vector<SpecialEdge>& sp_edges,
    std::vector<SHEET_NODE>& sheet_nodes_sorted, bool is_debug = false) {
  if (is_debug) {
    logger().debug(
        "calling load_sheet_nodes_sorted_for_sp_corner for corner tag {} ...",
        sp_corner.mat_tag);
    // sp_corner.print_info();
  }
  const auto& adj_edges_sorted = sp_corner.adj_edges_sorted;
  int num_adj_edges = adj_edges_sorted.size();
  // loop sharp edges first
  // sheet_nodes_sorted size == sharp edges size
  for (int cur_idx = 0; cur_idx < num_adj_edges; cur_idx++) {
    auto& cur_spe = sp_edges.at(adj_edges_sorted.at(cur_idx));
    if (cur_spe.type == EdgeType::CCE) continue;
    // create leaf nodes with sheet_seg of 2 tangent planes
    int id = sheet_nodes_sorted.size();
    cur_spe.sheet_node_id = id;
    // logger().debug("looping sharp edge cur_spe: {}, type {}, sheet_node_id:
    // {}", cur_spe.id, cur_spe.type, id);
    int related_mat_tag = cur_spe.mat_se_tags[0] == sp_corner.mat_tag
                              ? cur_spe.mat_se_tags[1]
                              : cur_spe.mat_se_tags[0];
    SHEET_NODE node(id, 1 /*type*/, related_mat_tag,
                    sp_corner.mat_tag /*corner_tag*/);
    node.sp_edge_id = cur_spe.id;
    SHEET_SEG sheet_seg;  // store 2 tangent planes
    for (int i = 0; i < 2; i++) {
      // tangent points needs to be close to corner
      // using sp_corner.min_corner_len
      if (sp_corner.min_corner_len == DBL_MAX)
        log_and_throw("ERROR: sp_corner.min_corner_len cannot be DBL_MAX");
      Vector3 tan_point =
          sp_corner.pos +
          sp_corner.min_corner_len *
              get_direction(sp_corner.pos, cur_spe.adj_ref_tan_points[i]);
      TangentPlane tan_pl(cur_spe.adj_ref_normals[i], tan_point);
      sheet_seg.tan_pls.push_back(tan_pl);
    }
    node.sheet_segs.push_back(sheet_seg);
    sheet_nodes_sorted.push_back(node);
  }

  // loop concave edges if any
  // update the tangent element of leaf nodes
  for (int cur_idx = 0; cur_idx < num_adj_edges; cur_idx++) {
    const auto& cur_spe = sp_edges.at(adj_edges_sorted.at(cur_idx));
    // logger().debug("looping concave edge cur_spe: {}, type {}", cur_spe.id,
    //                cur_spe.type);
    if (cur_spe.type != EdgeType::CCE) continue;
    // cur_spe is a concave line
    std::vector<int> neigh_ids;  // check prev and next
    neigh_ids.push_back((cur_idx - 1) % num_adj_edges);
    neigh_ids.push_back((cur_idx + 1) % num_adj_edges);
    for (const auto& neigh_idx : neigh_ids) {
      // neigh_spe can be sharp line or concave line
      const auto& neigh_spe = sp_edges.at(adj_edges_sorted.at(neigh_idx));
      // if neigh_spe is concave line, then sheet_node_id == -1
      if (neigh_spe.sheet_node_id == -1) continue;
      auto& node = sheet_nodes_sorted.at(neigh_spe.sheet_node_id);
      // if neighboring any sharp edge, then update sheet
      if (neigh_spe.type != EdgeType::SE) continue;

      // create new sheet_seg of 1 tangent plane and 1 concave line
      SHEET_SEG sheet_seg;
      sheet_seg.tan_cc_line_ids.insert(cur_spe.cc_line_id);
      const auto& cur_cc_line = tan_cc_lines.at(cur_spe.cc_line_id);
      for (int i = 0; i < 2; i++) {
        const auto& normal = neigh_spe.adj_ref_normals[i];
        if (cur_cc_line.is_normal_covered_by_adj_fs(normal, esp_degree_10))
          continue;
        TangentPlane tan_pl(normal, neigh_spe.adj_ref_tan_points[i]);
        sheet_seg.tan_pls.push_back(tan_pl);
      }
      node.sheet_segs.push_back(sheet_seg);

      // if neighboring sharp edge already has concave edges
      // then create new sheet_seg of 2 concave lines
      for (const auto& sheet_seg1 : node.sheet_segs) {
        if (sheet_seg1.tan_cc_line_ids.size() > 0 &&
            sheet_seg1.tan_cc_line_ids.find(cur_spe.cc_line_id) ==
                sheet_seg1.tan_cc_line_ids.end()) {
          SHEET_SEG sheet_seg2;
          sheet_seg2.tan_cc_line_ids.insert(sheet_seg1.tan_cc_line_ids.begin(),
                                            sheet_seg1.tan_cc_line_ids.end());
          sheet_seg2.tan_cc_line_ids.insert(cur_spe.cc_line_id);
          node.sheet_segs.push_back(sheet_seg2);
        }
      }  // for node.sheet_segs
    }    // for prev and next sp_edges
  }      // for num_adj_edges

  if (is_debug) {
    for (const auto& node : sheet_nodes_sorted) {
      node.print_info();
    }
  }
}

// union or diff tangent planes of sheet_seg1 + sheet_seg2
// and filter out those tangent planes coverd by tangent concave lines
// save the filtered tangent planes to new_sheet_seg
void union_or_diff_tan_pls_and_filter_by_tan_cc_lines(
    const bool is_union, const std::vector<TangentConcaveLine>& tan_cc_lines,
    const SHEET_SEG& sheet_seg1, const SHEET_SEG& sheet_seg2,
    SHEET_SEG& new_sheet_seg) {
  std::vector<TangentPlane> tmp_tan_pls;
  if (is_union)
    TangentPlane::union_two_tan_pls_vectors(sheet_seg1.tan_pls,
                                            sheet_seg2.tan_pls, tmp_tan_pls);
  else
    TangentPlane::diff_two_tan_pls_vectors(sheet_seg1.tan_pls,
                                           sheet_seg2.tan_pls, tmp_tan_pls);

  // filter tangent planes in tmp_tan_pls that are covered by tan_cc_lines
  // save to new_sheet_seg.tan_pls is not covered by cc_lines
  for (const auto& tan_pl : tmp_tan_pls) {
    bool is_tan_pl_covered = false;
    for (const auto& cc_line_id : new_sheet_seg.tan_cc_line_ids) {
      const auto& cur_cc_line = tan_cc_lines.at(cc_line_id);
      if (cur_cc_line.is_normal_covered_by_adj_fs(tan_pl.normal,
                                                  esp_degree_10)) {
        is_tan_pl_covered = true;
        continue;
      }
      if (is_tan_pl_covered) break;
    }
    if (!is_tan_pl_covered) new_sheet_seg.tan_pls.push_back(tan_pl);
  }  // for tmp_tan_pls
}

// 1. update seam_segs
// 2. update new_sheet_segs
// each seam_segs can grow a new sheet segment
// so seam_segs.size == new_sheet_segs.size
bool load_tangent_info_given_two_sheet_nodes(
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    const SHEET_NODE& node1, const SHEET_NODE& node2,
    std::vector<SHEET_SEG>& seam_segs, std::vector<SHEET_SEG>& new_sheet_segs,
    bool is_debug = false) {
  if (is_debug)
    logger().debug("calling load_tangent_info_given_two_sheet_nodes...");

  seam_segs.clear();
  new_sheet_segs.clear();
  for (const auto& sheet_seg1 : node1.sheet_segs) {
    for (const auto& sheet_seg2 : node2.sheet_segs) {
      // union tangent info of two sheet segments
      // for computing seam sphere
      // e.g. (sheet_seg1 || sheet_seg2)
      SHEET_SEG one_seam_seg;
      // diff tangent info of a new sheet segment
      // each seam sphere can grow a new sheet segment
      // e.g. (sheet_seg1 || sheet_seg2) - (sheet_seg1 && sheet_seg2)
      SHEET_SEG one_new_sheet_seg;

      // union and diff tangent concave lines if any
      one_seam_seg.tan_cc_line_ids = sheet_seg1.tan_cc_line_ids;
      one_seam_seg.tan_cc_line_ids.insert(sheet_seg2.tan_cc_line_ids.begin(),
                                          sheet_seg2.tan_cc_line_ids.end());
      set_difference(sheet_seg1.tan_cc_line_ids, sheet_seg2.tan_cc_line_ids,
                     one_new_sheet_seg.tan_cc_line_ids);

      // union and diff tangent planes if any
      // tan_cc_line_ids must given before calling
      union_or_diff_tan_pls_and_filter_by_tan_cc_lines(
          true /*union*/, tan_cc_lines, sheet_seg1, sheet_seg2, one_seam_seg);
      union_or_diff_tan_pls_and_filter_by_tan_cc_lines(
          false /*diff*/, tan_cc_lines, sheet_seg1, sheet_seg2,
          one_new_sheet_seg);

      if (is_debug) {
        logger().debug("[load_tangent_info] one_seam_seg:");
        one_seam_seg.print_info();
        logger().debug("[load_tangent_info] one_new_sheet_seg:");
        one_new_sheet_seg.print_info();
      }

      bool is_good = false;
      // sanity check of one_seam_seg
      int size_tan_cc_lines = one_seam_seg.tan_cc_line_ids.size();
      int size_tan_pls = one_seam_seg.tan_pls.size();
      if (size_tan_cc_lines + size_tan_pls >= 3) {
        is_good = true;
      }
      // sanity check of one_new_sheet_seg
      size_tan_cc_lines = one_new_sheet_seg.tan_cc_line_ids.size();
      size_tan_pls = one_new_sheet_seg.tan_pls.size();
      if (size_tan_cc_lines + size_tan_pls == 2) {
        is_good = is_good && true;
      }
      if (is_good) {
        seam_segs.push_back(one_seam_seg);
        new_sheet_segs.push_back(one_new_sheet_seg);
      }

      if (is_debug) {
        logger().debug("[load_tangent_info] is_good: {}", is_good);
      }
    }  // for node2.sheet_segs
  }    // for node1.sheet_segs

  if (is_debug)
    logger().debug(
        "[load_tangent_info] seam_segs size: {}, new_sheet_segs size: {}",
        seam_segs.size(), new_sheet_segs.size());

  return seam_segs.size() > 0 && new_sheet_segs.size() > 0 &&
         seam_segs.size() == new_sheet_segs.size();
}

// helper function
// for nodes left, we store added_for_two_spheres
// to sphere with is_added_for_corner == true (type == 2)
void post_process_sheet_tree_to_merge(
    std::queue<std::unique_ptr<SHEET_NODE>>& sheet_queue,
    std::vector<MVertex>& all_medial_spheres, bool is_debug = false) {
  if (sheet_queue.size() <= 1) return;  // do nothing
  if (is_debug)
    logger().debug(
        "calling post_process_sheet_tree_to_merge with sheet_queue size "
        "{}...",
        sheet_queue.size());

  std::unique_ptr<SHEET_NODE> node_mom = nullptr;
  std::vector<std::unique_ptr<SHEET_NODE>> node_children;

  while (!sheet_queue.empty()) {
    auto node = std::move(sheet_queue.front());  // not a copy
    sheet_queue.pop();
    // find first node with type == 2
    if (node->type == 2 && node_mom == nullptr)
      node_mom = std::move(node);
    else
      node_children.push_back(std::move(node));
  }
  if (node_mom == nullptr) {
    if (is_debug)
      logger().debug(
          "cannot find seam node with type 2, node_mom is nullptr, do noting "
          "for this corner");
    return;
  }

  // for each node_child, save mat face as (node_mon, corner, node_child)
  int corner_tag = node_mom->corner_tag;
  int node_mom_tag = node_mom->related_mat_tag;
  auto& mat_p_mon = all_medial_spheres.at(node_mom_tag);
  for (const auto& node_child : node_children) {
    int node_child_tag = node_child->related_mat_tag;
    auto& mat_p_child = all_medial_spheres.at(node_child_tag);
    mat_p_mon.added_for_two_spheres.insert({{corner_tag, mat_p_child.tag}});
  }
  if (is_debug)
    logger().debug("create mat face {} with each pair {}", mat_p_mon.tag,
                   mat_p_mon.added_for_two_spheres);
}

// TODO: change new_mat_p.tan_cc_lines to set of ids
void grow_sheet_tree_for_one_sp_corner(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper, const SpecialCorner& sp_corner,
    const std::vector<TangentConcaveLine> tan_cc_lines,
    const std::vector<SpecialEdge>& sp_edges,
    std::vector<MVertex>& all_medial_spheres,
    std::vector<SHEET_NODE>& sheet_nodes_sorted,
    std::queue<std::unique_ptr<SHEET_NODE>>& sheet_queue,
    bool is_debug = false) {
  if (is_debug)
    logger().debug(
        "[sp_corner] start growing sheet tree for sp_corner: {} with tag {}",
        sp_corner.id, sp_corner.mat_tag);

  // load sorted sheet nodes to a queue
  for (const SHEET_NODE& sheet : sheet_nodes_sorted) {
    sheet_queue.push(std::make_unique<SHEET_NODE>(sheet));
  }

  // loop sheet_queue
  const int num_while_loop = 20;
  int num_loop = 0;
  bool is_tree_good = false;
  std::vector<SHEET_SEG> seam_segs, new_sheet_segs;
  // stop when only two nodes left
  while (sheet_queue.size() > 2 && num_loop < num_while_loop) {
    if (is_debug)
      logger().debug("[sp_corner] looping: sheet_queue size {}, num_loop: {}",
                     sheet_queue.size(), num_loop);
    num_loop++;
    // Also important is that after the move, the top element in the queue is
    // a unique_ptr equal to nullptr. The pop is needed to remove this "empty"
    // unique_ptr from the queue
    auto cur_node = std::move(sheet_queue.front());
    sheet_queue.pop();
    const auto& next_node = sheet_queue.front();  // not a copy

    if (is_debug) {
      logger().debug("[sp_corner] cur_node: {}, next_node: {}", cur_node->id,
                     next_node->id);
      cur_node->print_info();
      next_node->print_info();
    }
    seam_segs.clear();
    new_sheet_segs.clear();
    is_tree_good = load_tangent_info_given_two_sheet_nodes(
        tan_cc_lines, *cur_node, *next_node, seam_segs, new_sheet_segs,
        is_debug /*is_debug*/);

    // only 2 nodes left, this is fine
    if (!is_tree_good) {
      sheet_queue.push(std::move(cur_node));
      break;  // break while queue
    }

    // for each seam seg, try to create seam sphere
    // once created, we do not run the rest seam segments
    bool is_seam_sphere_created = false;
    MVertex new_mat_p(-1, 0., 0., 0., DBL_MAX);
    for (const auto& one_seam_seg : seam_segs) {
      if (is_debug) {
        logger().debug("---- creating new seam sphere ...");
        logger().debug(
            "[sp_corner one_seam_seg] has tan_pls {}, tan_cc_lines {}",
            one_seam_seg.tan_pls.size(), one_seam_seg.tan_cc_line_ids.size());
        // one_seam_seg.print_info();
      }

      new_mat_p.tan_cc_lines.clear();
      for (const int& cc_line_id : one_seam_seg.tan_cc_line_ids) {
        new_mat_p.tan_cc_lines.push_back(tan_cc_lines.at(cc_line_id));
      }
      new_mat_p.tan_planes = one_seam_seg.tan_pls;
      // force to restrict to tangent point, so alpha1 = 0.5 now
      double alpha1 = 0.5, alpha2 = 1., alpha3 = 1.;
      double break_threshold = SCALAR_ZERO_N3;
      int iteration_limit = 20;  // iterate less
      bool is_success = matfp::iterate_sphere(
          sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, new_mat_p,
          is_debug /*is_debug*/, alpha1, alpha2, alpha3, break_threshold,
          iteration_limit);
      if (!is_success) {
        if (is_debug) {
          logger().debug(
              "[sp_corner one_seam_seg] cannot create new seam sphere");
          new_mat_p.print_tan_planes();
          new_mat_p.print_cc_lines();
        }
        continue;
      }

      // found intersection, store connectivity
      new_mat_p.tag = all_medial_spheres.size();
      new_mat_p.is_added_for_corner = true;
      new_mat_p.is_sphere_updated = true;
      new_mat_p.added_for_two_spheres.insert(
          {{cur_node->related_mat_tag, cur_node->corner_tag}});
      new_mat_p.added_for_two_spheres.insert(
          {{next_node->related_mat_tag, next_node->corner_tag}});
      is_seam_sphere_created = true;
      all_medial_spheres.push_back(new_mat_p);
      // grow a new sheet node
      // update tangent info of new_sheet_node
      SHEET_NODE new_sheet_node(sheet_nodes_sorted.size(), 2 /*type*/,
                                new_mat_p.tag /*related_mat_tag*/,
                                sp_corner.mat_tag /*corner_tag*/);
      // new_sheet_node has sheet_segs of all possible combinations
      // not just one segment from new_sheet_segs that matching one_seam_seg
      // (store all possible combinations, not new_sheet_segs[i])
      new_sheet_node.sheet_segs = new_sheet_segs;
      // new_sheet_node.children.push_back(cur_node);
      // new_sheet_node.children.push_back(next_node);
      sheet_nodes_sorted.push_back(new_sheet_node);
      sheet_queue.pop();  // pop next_node as well
      sheet_queue.push(std::make_unique<SHEET_NODE>(sheet_nodes_sorted.back()));

      if (is_debug) {
        logger().debug(
            "[sp_corner one_seam_seg] sp_corner.id {} with "
            "sp_corner.mat_tag: "
            "{} created "
            "new seam sphere {} of sheet node {}",
            sp_corner.id, sp_corner.mat_tag, new_mat_p.tag, new_sheet_node.id);
        new_mat_p.print_info();
      }
      break;
    }  // for seam_segs

    // if cur_node and next_node CANNOT find intersection
    if (!is_seam_sphere_created) {
      // move to the back of the queue
      sheet_queue.push(std::move(cur_node));
    }
  }  // while sheet_queue

  if (!is_tree_good) {
    logger().debug("[sp_corner] ERROR: sheet_queue is not good with size {}",
                   sheet_queue.size());
    log_and_throw("ERROR");
  }

  // update mat_p.added_for_two_spheres for rest of nodes
  post_process_sheet_tree_to_merge(sheet_queue, all_medial_spheres, is_debug);

  if (is_debug)
    logger().debug("[sp_corner] sheet_queue {}, is_tree_good: {}",
                   sheet_queue.size(), is_tree_good);
}

////////////////////////////////////////////////////////////////////////////////
//  Main Functions
////////////////////////////////////////////////////////////////////////////////
//
// 1. remove spheres too close to sharp edges
// 2. add new sphere for preserving corners
void prepare_all_medial_spheres(ThreeDimensionalShape& shape3D, bool is_debug) {
  logger().debug("preparing all_medial_spheres ...");
  //////////
  // II: Weighted Delaunay Triangulation
  // 4 dimensional
  int dim = 4;
  std::vector<MVertex>& all_medial_spheres = shape3D.all_medial_spheres;
  std::vector<int>& valid_medial_spheres = shape3D.valid_medial_spheres;
  std::map<int, int>& all_to_valid_medial_spheres =
      shape3D.all_to_valid_medial_spheres;

  std::vector<SpecialEdge>& sp_edges = shape3D.sp_edges;
  std::vector<SpecialCorner>& sp_corners = shape3D.sp_corners;

  const GEO::Mesh& sf_mesh = shape3D.sf_mesh;
  const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross =
      shape3D.ref_fs_pairs_not_cross;
  const AABBWrapper& aabb_wrapper = shape3D.aabb_wrapper;
  const std::vector<TangentConcaveLine>& tan_cc_lines = shape3D.tan_cc_lines;

  ///////////
  // Sharp Edge: Deletion
  // TODO: use kd_tree
  remove_medial_spheres_too_close_to_sharp_edges(
      shape3D.se_spheres, all_medial_spheres, true /*is_using_dilated_radius*/,
      is_debug /*is_debug*/);

  //////////
  // Corners
  // add non-feature MAT point close to corners but not too close
  // we lock these newly added mat points during optimization
  for (int i = 0; i < sp_corners.size(); i++) {
    // if (i != 3) continue;
    auto& one_sp_corner = sp_corners[i];
    // update corner small region first
    load_corners_min_length(sp_edges, one_sp_corner, all_medial_spheres);

    std::vector<SHEET_NODE> sheet_nodes_sorted;
    load_sheet_nodes_sorted_for_sp_corner(one_sp_corner, tan_cc_lines, sp_edges,
                                          sheet_nodes_sorted, is_debug);
    // sheet_queue should be size 1 after growing
    // we only store the root node for no use
    // everything we need has been stored in MVertex::added_for_two_spheres
    std::queue<std::unique_ptr<SHEET_NODE>> sheet_queue;
    grow_sheet_tree_for_one_sp_corner(
        sf_mesh, ref_fs_pairs_not_cross, aabb_wrapper, one_sp_corner,
        tan_cc_lines, sp_edges, all_medial_spheres, sheet_nodes_sorted,
        sheet_queue, is_debug);
  }
}

// For each given non-feature sphere, we check K = 5 nearest sharp edges
// and add new feature spheres as much as possible.
// If not, then delete this non-feature sphere.
bool add_or_delete_for_se(std::vector<MVertex>& all_medial_spheres,
                          std::set<std::array<int, 2>>& se_spheres,
                          std::map<int, int>& se_kd_tree_idx_to_se_tag,
                          GEO::NearestNeighborSearch_var& se_kd_tree,
                          std::vector<double>& se_kd_points,
                          bool is_check_updated_only,
                          bool is_using_dilated_radius, bool is_debug) {
  if (is_debug)
    logger().debug(
        "[se_rec] calling add_or_delete_for_se with se_spheres: {}, "
        "is_debug: "
        "{}...",
        se_spheres.size(), is_debug);
  std::vector<MVertex> new_feature_spheres, all_new_feature_spheres;
  std::set<std::array<int, 2>> new_se_spheres;
  std::set<int> deleted_spheres, new_spheres_tags;
  std::map<int, std::unordered_set<int>> se_spheres_adjs;

  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) continue;
    if (mat_p.is_a_feature_sphere()) continue;
    // only check updated spheres
    if (is_check_updated_only && !mat_p.is_sphere_updated) {
      continue;
    }
    // reset mat_p.is_sphere_updated
    mat_p.is_sphere_updated = false;
    new_feature_spheres.clear();
    std::set<std::array<int, 2>> se_spheres_copy = se_spheres;
    int latest_tag = all_medial_spheres.size();
    double sq_radius = mat_p.sq_radius;
    if (is_using_dilated_radius) {
      /* use dilated radius !!!!*/
      mat_p.dilate_sphere_radius();
      sq_radius = mat_p.sq_radius_dilated;
    }
    // fetch K nearest neighbors
    const size_t K_nb_neigh = 5;
    std::set<std::array<int, 2>> k_nearest_se;
    if (se_spheres_adjs.empty()) {
      reload_se_spheres_adjs(se_spheres, se_spheres_adjs);
    }
    if (se_kd_tree == nullptr ||
        se_kd_tree->nb_points() != se_spheres_adjs.size()) {
      reload_se_kd_tree(all_medial_spheres, se_spheres, se_spheres_adjs,
                        se_kd_tree, se_kd_points, se_kd_tree_idx_to_se_tag,
                        is_debug);
    }
    fetch_k_nearest_sharp_edges(K_nb_neigh, all_medial_spheres, se_spheres_adjs,
                                se_kd_tree, se_kd_tree_idx_to_se_tag, mat_p.pos,
                                k_nearest_se, is_debug);
    if (is_debug)
      logger().debug("[se_rec] mat_p {} check sharp edges: {}", mat_p.tag,
                     k_nearest_se);
    for (const auto one_se : k_nearest_se) {
      auto& se_mat1 = all_medial_spheres.at(one_se[0]);
      auto& se_mat2 = all_medial_spheres.at(one_se[1]);

      // if (mat_p.tag == 43432 &&
      //     (se_mat1.tag == 36324 || se_mat2.tag == 36324)) {
      //   is_debug = true;
      //   logger().debug("mat {} checking sharp edge ({},{})", mat_p.tag,
      //                  se_mat1.tag, se_mat2.tag);
      // }
      ///////////////////
      if (!is_rpd_break_one_sharp_edge(se_mat1.pos, se_mat2.pos, mat_p.pos,
                                       sq_radius, is_debug /*is_debug*/)) {
        // is_debug = false;
        continue;
      }
      // is_debug = false;

      ///////////////////
      // reaching the small region of corners
      // since we do not want to change
      // anything in that region, so we delete the sphere once break
      if (se_mat1.is_on_corner || se_mat2.is_on_corner) {
        // if ((se_mat1.is_on_corner && se_mat1.min_corner_len == DBL_MAX) ||
        //     (se_mat2.is_on_corner && se_mat2.min_corner_len == DBL_MAX))
        //   log_and_throw(
        //       "ERROR: se_mat1 or se_mat2 min_corner_len cannot be DBL_MAX");
        if (mat_p.is_added_for_corner) continue;  // do not delete
        // double dist = get_distance_between_two_vectors(mat_p.pos,
        // se_mat1.pos); double dist_threshold = se_mat1.min_corner_len; if
        // (se_mat2.is_on_corner) {
        //   dist = get_distance_between_two_vectors(mat_p.pos, se_mat2.pos);
        //   dist_threshold = se_mat2.min_corner_len;
        // }
        // if (dist <= dist_threshold * 2) {
        mat_p.is_deleted = true;
        deleted_spheres.insert(mat_p.tag);
        // if (is_debug)
        //   logger().debug(
        //       "[se_rec] sphere {} reaches and breaks corner region with "
        //       "dist {}, dist_threshold: {}, delete",
        //       mat_p.tag, dist, dist_threshold);
        continue;
        // }
      }

      ///////////////////
      // check if adding new feature spheres
      if (is_debug)
        logger().debug(
            "[se_rec] mat_p {} breaks se ({},{}) going to add more feature "
            "spheres",
            mat_p.tag, se_mat1.tag, se_mat2.tag);

      ///////////
      // latest_tag will also be updated
      int latest_tag_prev = latest_tag;
      new_feature_spheres.clear();
      new_se_spheres.clear();
      bool is_sucess = check_and_add_to_se_recursively(
          {{se_mat1, se_mat2}}, mat_p, new_feature_spheres, new_se_spheres,
          latest_tag, is_using_dilated_radius, false /*is_debug*/);
      // adding new feature spheres cannot help
      // delete and reset
      if (!is_sucess) {
        if (is_debug)
          logger().debug("[se_rec] cannot add spheres for mat {}, delete",
                         mat_p.tag);
        mat_p.is_deleted = true;
        latest_tag = latest_tag_prev;
        deleted_spheres.insert(mat_p.tag);
        continue;
      }
      // nothing change
      if (new_se_spheres.size() == 1 || new_feature_spheres.empty()) {
        continue;
      }

      /////////
      // erase se_spheres and mat.se_adj_se
      se_spheres.erase(one_se);
      int num_rm1 = se_mat1.se_adj_se.erase(se_mat2.tag);
      int num_rm2 = se_mat2.se_adj_se.erase(se_mat1.tag);
      if (is_debug && (num_rm1 == 0 || num_rm2 == 0)) {
        logger().error("[se_rec] se_mat1 {} num_rm1 {}, se_mat2 {} num_rm2 {}",
                       se_mat1.tag, num_rm1, se_mat2.tag, num_rm2);
      }
      /////////
      // update new_se_spheres and all_new_feature_spheres
      se_spheres.insert(new_se_spheres.begin(), new_se_spheres.end());
      for (const auto& new_mat : new_feature_spheres) {
        mat_p.newly_added_feature_mat_tags.insert(new_mat.tag);
        all_new_feature_spheres.push_back(new_mat);
        new_spheres_tags.insert(new_mat.tag);
      }
    }  // for k_nearest_se

    if (is_debug && all_new_feature_spheres.size() > 0)
      logger().debug("[se_rec] mat_p {} newly inserted {} feature spheres",
                     mat_p.tag, all_new_feature_spheres.size());

    // sanity check
    if (all_medial_spheres.size() + all_new_feature_spheres.size() !=
        latest_tag) {
      logger().error(
          "[se_rec] ERROR, mat_: {} has all_medial_spheres {} + "
          "all_new_feature_spheres {} "
          "!= "
          "latest_tag {}",
          mat_p.tag, all_medial_spheres.size(), all_new_feature_spheres.size(),
          latest_tag);
      log_and_throw("ERROR");
    }
    // save all_new_feature_spheres
    for (const auto& new_mat : all_new_feature_spheres) {
      all_medial_spheres.push_back(new_mat);
    }
    if (all_new_feature_spheres.size() > 0) {
      // reload se_kd_tree
      reload_se_spheres_adjs(se_spheres, se_spheres_adjs);
      reload_se_kd_tree(all_medial_spheres, se_spheres, se_spheres_adjs,
                        se_kd_tree, se_kd_points, se_kd_tree_idx_to_se_tag,
                        is_debug);
    }
    all_new_feature_spheres.clear();
  }  // for all_medial_spheres

  // update mat.se_adj_se by se_spheres
  for (const auto& se : se_spheres) {
    all_medial_spheres[se[0]].se_adj_se.insert(se[1]);
    all_medial_spheres[se[1]].se_adj_se.insert(se[0]);
  }

  if (is_debug) {
    logger().debug("[se_rec] new_spheres_tags {}: {}", new_spheres_tags.size(),
                   new_spheres_tags);
    logger().debug("[se_rec] deleted_spheres {}: {}", deleted_spheres.size(),
                   deleted_spheres);
  }
}

}  // namespace matfp