// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/MeshFeatureProcessor.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////
// Helper Functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Sharp Edges
int push_new_sample(std::vector<double>& points, const Vector3& pm,
                    int dim = 3) {
  int seed_idx = points.size() / dim;  // default dimension = 3
  points.push_back(pm[0]);
  points.push_back(pm[1]);
  points.push_back(pm[2]);
  if (dim == 4) points.push_back(0);  // add 0 as radius
  return seed_idx;
}

void add_new_sample_to_sharp_edges_and_lock(
    const Vector3& p0, const Vector3& p1, const int idx0, const int idx1,
    const double ideal_length, std::vector<double>& points,
    std::set<int>& point_is_locked,
    std::vector<std::array<int, 2>>& newly_added_se, int dim = 3) {
  newly_added_se.clear();
  Vector3 vector = p1 - p0;
  int num_new_points = std::ceil(vector.norm() / ideal_length) - 1;
  // logger().debug("going to add {} new points to edge ({},{}) with vector len
  // {}, ideal len {}",
  //     num_new_points, idx0, idx1,
  //     vector.norm(), ideal_length
  // );
  std::array<int, 2> e;
  if (num_new_points <= 0) {
    // logger().debug("No new points should be added to sharp edge ({},{})",
    // idx0, idx1); add original edge, but {{idx0, idx1}} are mapped to indices
    // of seed_points
    e = {{idx0, idx1}};
    std::sort(e.begin(), e.end());
    newly_added_se.push_back(e);
    return;
  }
  Vector3 step = vector / (num_new_points + 1);
  Vector3 pm = p0;
  int new_seed_idx = -1;
  for (int i = 0; i < num_new_points; i++) {
    pm += step;
    new_seed_idx = push_new_sample(points, pm, dim);
    point_is_locked.insert(new_seed_idx);
    // update newly_added_se
    if (i == 0) {
      e = {{idx0, new_seed_idx}};
    } else {
      e = {{new_seed_idx - 1, new_seed_idx}};
    }
    std::sort(e.begin(), e.end());
    newly_added_se.push_back(e);
  }
  // add last pair
  if (new_seed_idx != -1) {
    e = {{new_seed_idx, idx1}};
    std::sort(e.begin(), e.end());
    newly_added_se.push_back(e);
  }
}

void store_special_edge(const MeshWrapper& sf_mesh_wrapper,
                        const std::vector<double>& seed_points,
                        const std::array<int, 2>& old_e,
                        const std::array<int, 2>& new_e, const EdgeType& type,
                        std::vector<SpecialEdge>& sp_edges,
                        std::map<std::array<int, 2>, int>& map_to_sp_edges) {
  int sp_eid = sp_edges.size();
  map_to_sp_edges[new_e] = sp_eid;

  SpecialEdge sp_e(sp_eid, type, new_e);
  Vector3 p0(seed_points[new_e[0] * 3], seed_points[new_e[0] * 3 + 1],
             seed_points[new_e[0] * 3 + 2]);
  Vector3 p1(seed_points[new_e[1] * 3], seed_points[new_e[1] * 3 + 1],
             seed_points[new_e[1] * 3 + 2]);
  sp_e.vs_pos = {{p0, p1}};
  std::vector<int> fs_pair;
  std::vector<Vector3> fs_normals;
  // old_e is mapping to sf_mesh
  // new_e is mapping to seed_points!!!!
  sf_mesh_wrapper.get_fs_pair_given_edge(old_e, fs_pair);
  fs_normals.push_back(sf_mesh_wrapper.get_f_normal(fs_pair[0]));
  fs_normals.push_back(sf_mesh_wrapper.get_f_normal(fs_pair[1]));
  sp_e.adj_ref_fs_pair = {{fs_pair[0], fs_pair[1]}};
  sp_e.adj_ref_normals = {{fs_normals[0], fs_normals[1]}};
  sp_e.adj_ref_tan_points = {{sf_mesh_wrapper.get_f_centroid(fs_pair[0]),
                              sf_mesh_wrapper.get_f_centroid(fs_pair[1])}};
  sp_edges.push_back(sp_e);
}

////////////////////////////////////////////////////////////////////////////////
// Corners
Vector3 get_c_normal(const SpecialCorner& sp_corner,
                     const std::vector<SpecialEdge>& sp_edges) {
  Vector3 c_normal(0, 0, 0);
  int n_count = 0;  // init to 0 is important!! (don't ask me why)
  for (const int sp_eid : sp_corner.adj_edges) {
    const SpecialEdge& sp_e = sp_edges.at(sp_eid);
    c_normal += sp_e.adj_ref_normals[0];
    c_normal += sp_e.adj_ref_normals[1];
    n_count += 2;
  }
  c_normal /= n_count;
  return c_normal;
  ///////////
  // If global orientation of triangles are counter-clockwise
  // then we do not need to do anything. But if clockwise,
  // multiply -1 can give us proper normal pointing outside.
  // not that -1 here is just for pretty drawing, not contribute to the
  // calculation c_normal *= -1.;
  //
  // logger().info("corner {} has normal ({},{},{}), n_count {}",
  //     pair.first, c_normal[0], c_normal[1], c_normal[2], n_count);
}

void sort_sp_corner_adj_sp_edges(const std::vector<SpecialEdge>& sp_edges,
                                 SpecialCorner& one_sp_corner) {
  const std::set<int>& c_adj_edges = one_sp_corner.adj_edges;
  if (c_adj_edges.empty() || c_adj_edges.empty()) return;
  std::vector<int> c_adj_edges_sorted;
  std::vector<Vector3> c_adj_neighbors_pos;
  for (const int& adj_eid : c_adj_edges) {
    c_adj_edges_sorted.push_back(adj_eid);  // init
    const auto& sp_edge = sp_edges.at(adj_eid);
    for (int i = 0; i < 2; i++) {
      int seed_idx = sp_edge.vs_seed_ids[i];
      if (seed_idx == one_sp_corner.seed_id) continue;
      c_adj_neighbors_pos.push_back(sp_edge.vs_pos[i]);
    }
  }  // for c_adj_edges
  int num_neighs = c_adj_edges_sorted.size();
  if (num_neighs != c_adj_neighbors_pos.size()) {
    logger().error("num_neighs: {} != c_adj_neighbors_pos: {}", num_neighs,
                   c_adj_neighbors_pos.size());
    log_and_throw("ERROR");
  }

  // For each neighbor of corner, we find its projection on plane P.
  // Plane P is defined by point corner and normal corner_normal
  const Vector3 corner_pos = one_sp_corner.pos;
  const Vector3 corner_normal = get_c_normal(one_sp_corner, sp_edges);
  // logger().debug("corner {} has pos ({},{},{}) and c_normal ({},{},{})",
  //                one_sp_corner.id, corner_pos[0], corner_pos[1],
  //                corner_pos[2], corner_normal[0], corner_normal[1],
  //                corner_normal[2]);

  std::vector<std::pair<double, int>> angle_idx_pair;
  int neigh0_idx = c_adj_edges_sorted[0];  // first
  Vector3 neigh0_proj(0, 0, 0);
  double _1;
  for (int i = 0; i < num_neighs; i++) {
    if (i == 0) {
      matfp::project_point_on_plane(one_sp_corner.pos, corner_normal,
                                    c_adj_neighbors_pos[i], neigh0_proj, _1);
    } else {
      // Compute angle between two vectors embedded in a 3D plane
      // https://stackoverflow.com/a/16544330
      Vector3 neigh1_proj(0, 0, 0);
      matfp::project_point_on_plane(corner_pos, corner_normal,
                                    c_adj_neighbors_pos[i], neigh1_proj, _1);
      Vector3 vec0 = neigh0_proj - corner_pos;
      Vector3 vec1 = neigh1_proj - corner_pos;
      double dot = vec0.dot(vec1);                       // dx
      double det = corner_normal.dot(vec0.cross(vec1));  // dy
      // Since atan2 is discontinuity on -Pi and Pi
      // the counter clockwise order should be (0 to Pi) then (-Pi to 0)
      // therefore we move range (-Pi to 0) to (Pi to 2Pi) by adding 2Pi
      // then the order would be (0 to 2Pi)
      double angle = std::atan2(det, dot);
      if (angle < 0) {
        angle += 2 * PI;
      }
      // double angle = pseudoangle(det, dot);
      angle_idx_pair.push_back(std::make_pair(angle, c_adj_edges_sorted[i]));
    }
  }

  // Sorting the vector elements on the basis of
  // FIRST element (angle) of pairs in ascending/descending order.
  // (together with shape3D.corners_normals defined in find_edge_features())
  //
  // When global triangle (corner, neigh0, neigh1) is counter-clockwise,
  // sorting in ascending order would make sure
  // projected triangle (corner, neigh0, neigh1) are counter-clockwise,
  // vice versa.
  //
  // logger().debug("angle_idx_pair before sorting: {} ", angle_idx_pair);
  std::sort(angle_idx_pair.begin(), angle_idx_pair.end());
  // std::sort(angle_idx_pair.begin(), angle_idx_pair.end(), sort_in_rev);
  // logger().debug("angle_idx_pair after sorting: {} ", angle_idx_pair);

  if (angle_idx_pair.size() != num_neighs - 1)
    log_and_throw("ERROR: we ignore the first point on purpose");

  // reload
  c_adj_edges_sorted.clear();
  c_adj_edges_sorted.push_back(neigh0_idx);
  for (auto pair : angle_idx_pair) {
    c_adj_edges_sorted.push_back(pair.second);
  }
  one_sp_corner.adj_edges_sorted = c_adj_edges_sorted;
}

void update_corners_neighbors(const std::vector<SpecialEdge>& sp_edges,
                              std::vector<SpecialCorner>& sp_corners) {
  if (sp_edges.empty() || sp_corners.empty()) {
    logger().info(
        "ATTENTION: special edges or corners are zero or not be initialized");
    return;
  }
  // logger().debug("calling update_corners_neighbors ...");
  // re-sort
  for (auto& one_sp_corner : sp_corners) {
    sort_sp_corner_adj_sp_edges(sp_edges, one_sp_corner);
    // one_sp_corner.print_info();
  }
  // logger().debug("done update_corners_neighbors");
}

////////////////////////////////////////////////////////////////////////////////
// Main Functions
////////////////////////////////////////////////////////////////////////////////
//
// Sharp Edges
void load_sample_and_lock_sharp_edges(
    const MeshWrapper& sf_mesh_wrapper, const double& ideal_length,
    const double& s_edge_ideal_eps_rel,
    const GEO::NearestNeighborSearch_var& lfs_kd_tree,
    const std::vector<double>& lfs_min_values, const double r_sample,
    const std::set<int>& corners_set, std::vector<double>& seed_points,
    std::map<int, int>& seeds_map, std::set<int>& point_is_locked,
    std::vector<std::array<int, 2>>& s_edges,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals,
    std::vector<SpecialEdge>& sp_edges,
    std::map<std::array<int, 2>, int>& map_to_sp_edges, bool is_debug) {
  if (is_debug) logger().info("start load_sample_and_lock_sharp_edges...");
  std::vector<std::array<int, 2>> s_edges_new;
  std::map<std::array<int, 2>, std::array<Vector3, 2>> se_normals_new;
  // indices of new sharp edges by seed_points
  // -> original sharp edges by sf_mesh(sf_mesh_wrapper)
  std::map<std::array<int, 2>, std::array<int, 2>> se_new_to_old;
  //////////
  // Sharp edges
  // insert and lock feature seeds
  // edge {e[0], e[1]} is sorted based on their indices
  std::vector<std::array<int, 2>> newly_added_se;
  // will be updated by r-sample if defined LFS
  double se_ideal_length = ideal_length / 1.2;
  for (const auto old_se : s_edges) {
    // sanity check
    if (se_normals.find(old_se) == se_normals.end()) {
      logger().error("sharp edge {} cannot find normals", old_se);
      log_and_throw("ERROR");
    }
    // update se_ideal_length by LFS
    double minSizingField = DBL_MAX;
    for (int i = 0; i < 2; i++) {
      const Vector3& p_sf = sf_mesh_wrapper.get_v_pos(old_se[i]);
      std::vector<double> p = {{p_sf[0], p_sf[1], p_sf[2]}};
      GEO::index_t K_nb_neigh = 1;  // k nearest neighbors
      std::vector<GEO::index_t> k_neighbors(K_nb_neigh);
      std::vector<double> k_neighbors_sq_dist(K_nb_neigh);
      // LFS
      lfs_kd_tree->get_nearest_neighbors(
          K_nb_neigh, p.data(), k_neighbors.data(), k_neighbors_sq_dist.data());
      for (const GEO::index_t kn_id : k_neighbors) {
        double tmp = r_sample * lfs_min_values[kn_id];
        minSizingField = std::min(minSizingField, tmp);
      }
    }
    se_ideal_length = minSizingField * s_edge_ideal_eps_rel;
    // logger().debug("sharp edge {} has se_ideal_length: {}, minSizingField
    // {}", e, se_ideal_length, minSizingField);

    // // // // for debug, use constant
    // se_ideal_length = 100;

    //// lock two edge points if not corners
    int idx0 = old_se[0];
    Vector3 p0 = sf_mesh_wrapper.get_v_pos(idx0);
    auto it = seeds_map.find(idx0);
    if (it != seeds_map.end())
      idx0 = it->second;
    else
      idx0 = push_new_sample(seed_points, p0);
    point_is_locked.insert(idx0);
    seeds_map[old_se[0]] = idx0;

    int idx1 = old_se[1];
    Vector3 p1 = sf_mesh_wrapper.get_v_pos(idx1);
    auto it1 = seeds_map.find(idx1);
    if (it1 != seeds_map.end())
      idx1 = it1->second;
    else
      idx1 = push_new_sample(seed_points, p1);
    point_is_locked.insert(idx1);
    seeds_map[old_se[1]] = idx1;

    // if both sharp edge are corners, update se_ideal_length
    // make sure we have at least 1 sample in between
    if (corners_set.find(old_se[0]) != corners_set.end() &&
        corners_set.find(old_se[1]) != corners_set.end()) {
      double len = get_distance_between_two_vectors(p0, p1);
      se_ideal_length = min(se_ideal_length, len / 3.);
      if (is_debug)
        logger().debug(
            "sharp edge {} both on corners, update se_ideal_length {}", old_se,
            se_ideal_length);
    }

    // add new samples based on sizing field
    add_new_sample_to_sharp_edges_and_lock(p0, p1, idx0, idx1, se_ideal_length,
                                           seed_points, point_is_locked,
                                           newly_added_se);
    // updating se_normals and se_new_to_old
    std::array<Vector3, 2>& old_e_normals = se_normals[old_se];
    for (auto const new_se_sorted : newly_added_se) {
      se_normals_new[new_se_sorted] = old_e_normals;
      // mapping from new se to old se
      se_new_to_old[new_se_sorted] = old_se;
    }
    std::copy(newly_added_se.begin(), newly_added_se.end(),
              std::back_inserter(s_edges_new));
  }  // for s_edges done

  if (!s_edges_new.empty()) {
    if (is_debug)
      logger().debug("s_edges size: {}->{}", s_edges.size(),
                     s_edges_new.size());
    s_edges.clear();
    s_edges = s_edges_new;
  }
  if (!se_normals_new.empty()) {
    if (is_debug)
      logger().debug("se_normals: {}->{}", se_normals.size(),
                     se_normals_new.size());
    se_normals.clear();
    se_normals = se_normals_new;
  }

  // create special edges
  // now indices of s_edges are mapping to seed_points, not sf_mesh
  int old_sp_edges_size = sp_edges.size();
  for (const auto& new_se : s_edges) {
    const std::array<int, 2>& old_se = se_new_to_old.at(new_se);
    store_special_edge(sf_mesh_wrapper, seed_points, old_se, new_se,
                       EdgeType::SE, sp_edges, map_to_sp_edges);
  }
  if (is_debug)
    logger().debug("created {} special edges (SE):  {}->{}",
                   sp_edges.size() - old_sp_edges_size, old_sp_edges_size,
                   sp_edges.size());
  if (is_debug) logger().info("done load_sample_and_lock_sharp_edges");
}

// Concave Edges, will init tan_cc_lines
// TODO: also sample more points if not enough?
void load_concave_edges(const MeshWrapper& sf_mesh_wrapper,
                        std::vector<double>& seed_points,
                        std::map<int, int>& seeds_map,
                        std::vector<std::array<int, 2>>& cc_edges,
                        std::vector<SpecialEdge>& sp_edges,
                        std::vector<TangentConcaveLine>& tan_cc_lines,
                        std::map<std::array<int, 2>, int>& map_to_sp_edges,
                        bool is_debug) {
  // cc_edges is matching to sf_mesh now
  // but will be mathching to seed_points after this function
  if (is_debug) logger().info("start load_concave_edges...");
  std::vector<std::array<int, 2>> cc_edges_new;
  std::map<std::array<int, 2>, std::array<int, 2>> cce_new_to_old;
  for (const auto& cce_old : cc_edges) {
    int idx0 = cce_old[0];
    Vector3 p0 = sf_mesh_wrapper.get_v_pos(idx0);
    auto it = seeds_map.find(idx0);
    if (it != seeds_map.end())
      idx0 = it->second;
    else
      idx0 = push_new_sample(seed_points, p0);
    // do not lock
    seeds_map[cce_old[0]] = idx0;

    int idx1 = cce_old[1];
    Vector3 p1 = sf_mesh_wrapper.get_v_pos(idx1);
    auto it1 = seeds_map.find(idx1);
    if (it1 != seeds_map.end())
      idx1 = it1->second;
    else
      idx1 = push_new_sample(seed_points, p1);
    // do not lock
    seeds_map[cce_old[1]] = idx1;

    std::array<int, 2> cce_new = {{idx0, idx1}};
    std::sort(cce_new.begin(), cce_new.end());
    cc_edges_new.push_back(cce_new);
    cce_new_to_old[cce_new] = cce_old;
  }
  if (!cc_edges_new.empty()) {
    // logger().debug("cc_edges_new: {}", cc_edges_new);
    cc_edges.clear();
    cc_edges = cc_edges_new;
  }

  // create special edges and tan_cc_lines
  // now indices of s_edges are mapping to seed_points, not sf_mesh
  tan_cc_lines.clear();
  int old_sp_edges_size = sp_edges.size();
  for (int i = 0; i < cc_edges.size(); i++) {
    const std::array<int, 2>& new_cce = cc_edges[i];
    const std::array<int, 2>& old_cce = cce_new_to_old.at(new_cce);
    store_special_edge(sf_mesh_wrapper, seed_points, old_cce, new_cce,
                       EdgeType::CCE, sp_edges, map_to_sp_edges);

    int cc_line_id = tan_cc_lines.size();
    SpecialEdge& one_sp_edge = sp_edges.back();
    // create and store new TangentConcaveLine
    TangentConcaveLine new_cc_line(
        cc_line_id, {{one_sp_edge.vs_pos[0], one_sp_edge.vs_pos[1]}},
        one_sp_edge.adj_ref_normals, one_sp_edge.adj_ref_fs_pair,
        one_sp_edge.vs_seed_ids);
    tan_cc_lines.push_back(new_cc_line);

    // mapping sp_edges with tan_cc_lines
    one_sp_edge.cc_line_id = cc_line_id;
  }
  if (is_debug) {
    logger().debug(
        "created {} special edges (CCE), {}->{}, and tan_cc_lines {}",
        sp_edges.size() - old_sp_edges_size, old_sp_edges_size, sp_edges.size(),
        tan_cc_lines.size());
    logger().info("done load_concave_edges");
  }
}

// Corners Part1
void load_and_lock_corners(const MeshWrapper& sf_mesh_wrapper,
                           std::vector<double>& seed_points,
                           std::map<int, int>& seeds_map,
                           std::set<int>& point_is_locked,
                           std::vector<int>& corners, bool is_debug) {
  if (is_debug) logger().info("start load_and_lock_corners...");
  std::vector<int> corners_new;
  std::map<int, Vector3> corners_normals_new;
  ////////////
  // corners is matching sf_mesh now (old_c)
  // but will be matching to seed_points after this function (new_c)
  for (const auto& old_c : corners) {
    Vector3 p = sf_mesh_wrapper.get_v_pos(old_c);
    int new_c = -1;
    auto it = seeds_map.find(old_c);
    if (it != seeds_map.end())
      new_c = it->second;
    else
      new_c = push_new_sample(seed_points, p);
    point_is_locked.insert(new_c);
    seeds_map[old_c] = new_c;
    corners_new.push_back(new_c);
  }
  if (!corners_new.empty()) {
    corners.clear();
    corners = corners_new;
  }
  if (is_debug) logger().debug("updated corners: {}", corners.size());
  if (is_debug) logger().info("done load_and_lock_corners");
}

// Corners Part2
// NOTE:
// when calling, corners must mapping to seed_points, not sf_mesh!!!
void load_corner_special_edges(const std::vector<SpecialEdge>& sp_edges,
                               std::vector<int>& corners,
                               std::vector<SpecialCorner>& sp_corners,
                               bool is_debug) {
  if (is_debug) logger().debug("calling load_corner_special_edges...");
  sp_corners.clear();
  std::set<int> corners_set;
  corners_set.insert(corners.begin(), corners.end());
  std::map<int, int> map_seed_to_spc;  // seed_idx to SpecialCorner::id
  for (const auto& sp_e : sp_edges) {
    int sp_eid = sp_e.id;
    for (int i = 0; i < 2; i++) {
      // not a corner
      int seed_idx = sp_e.vs_seed_ids[i];
      Vector3 seed_pos = sp_e.vs_pos[i];
      if (corners_set.find(seed_idx) == corners_set.end()) continue;
      int spc_id = -1;
      if (map_seed_to_spc.find(seed_idx) != map_seed_to_spc.end()) {
        spc_id = map_seed_to_spc.at(seed_idx);
        SpecialCorner& sp_c = sp_corners[spc_id];
        sp_c.push_to_adj_edges(sp_eid);
      } else {
        // create new Corner
        spc_id = sp_corners.size();
        SpecialCorner sp_c(spc_id, seed_idx, seed_pos);
        map_seed_to_spc[seed_idx] = spc_id;
        sp_c.push_to_adj_edges(sp_eid);
        sp_corners.push_back(sp_c);
      }
    }  // for sp_e.vs_seed_ids
  }    // for sp_edges

  // update neighbors of each special corners
  update_corners_neighbors(sp_edges, sp_corners);
  if (is_debug) logger().debug("loaded special corners: {}", sp_corners.size());
  if (is_debug) logger().debug("done load_corner_special_edges");
}

}  // namespace matfp
