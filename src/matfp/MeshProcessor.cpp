// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/MeshProcessor.h"

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/points/kd_tree.h>
#include <geogram/points/nn_search.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>

#include "matfp/Common.h"
#include "matfp/MeshFeatureProcessor.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/external/Predicates.hpp"
#include "matfp/geogram/mesh/mesh_sampling.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////
// Bounding box
////////////////////////////////////////////////////////////////////////////////
void get_bb_corners(Parameters& params, const std::vector<double>& vertices,
                    Vector3& min, Vector3& max) {
  min = Vector3(vertices[0], vertices[1], vertices[2]);
  max = min;

  for (size_t j = 0; j < vertices.size(); j += 3) {
    for (int i = 0; i < 3; i++) {
      min(i) = std::min(min(i), vertices[j + i]);
      max(i) = std::max(max(i), vertices[j + i]);
    }
  }

  const Scalar dis = std::max(params.init_edge_length, params.eps_input * 2);
  for (int j = 0; j < 3; j++) {
    min[j] -= dis;
    max[j] += dis;
  }

  logger().debug("min = {} {} {}", min[0], min[1], min[2]);
  logger().debug("max = {} {} {}", max[0], max[1], max[2]);

  // params.bbox_min = min;
  // params.bbox_max = max;

  // params.bbox_min = Vector3(min[0], min[1], min[2]);
  // params.bbox_max = Vector3(max[0], max[1], max[2]);
}

////////////////////////////////////////////////////////////////////////////////
// Feature detection
////////////////////////////////////////////////////////////////////////////////
bool get_normals_given_edge(
    const std::array<int, 2>& e,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& edge_normals,
    std::array<Vector3, 2>& ns) {
  if (edge_normals.empty() || e.empty()) return false;

  // print_se_normals(edge_normals);

  // find normals
  bool is_found = false;
  for (int i = 0; i < 2; i++) {
    std::array<int, 2> tmp = {{e[i], e[(i + 1) % 2]}};
    // logger().debug("searching normals for edges {}", tmp);

    auto it = edge_normals.find(tmp);
    if (it != edge_normals.end()) {
      ns = it->second;
      is_found = true;
      break;
    }
  }
  if (!is_found) {
    logger().debug("Not found normals for edge {} (vertex order not matter)",
                   e);
    return false;
  }
  // logger().debug("Found edge {} has normal [({},{},{}), ({},{},{})]",
  // 	e,
  // 	ns[0][0], ns[0][1], ns[0][2],
  // 	ns[1][0], ns[1][1], ns[1][2]
  // );
  return true;
}

void print_se_normals(
    const std::map<std::array<int, 2>, std::array<Vector3, 2>>& edge_normals) {
  std::map<std::array<int, 2>, std::array<Vector3, 2>>::const_iterator it =
      edge_normals.begin();
  while (it != edge_normals.end()) {
    const std::array<int, 2> e = it->first;
    const std::array<Vector3, 2> ns = it->second;
    logger().info("Found edge {} has normal [({},{},{}), ({},{},{})]", e,
                  ns[0][0], ns[0][1], ns[0][2], ns[1][0], ns[1][1], ns[1][2]);
    it++;
  }
}

// sf_mesh is matching input_vertices & input_faces
void load_mesh_features_helper(
    const GEO::Mesh& sf_mesh, const std::vector<Vector3>& input_vertices,
    const std::vector<Vector3i>& input_faces,
    std::vector<std::array<int, 2>>& s_edges,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals,
    std::set<std::array<int, 2>>& se_ref_fs_pairs,
    std::vector<std::array<int, 2>>& cc_edges, std::vector<int>& corners,
    std::map<int, std::unordered_set<int>>& conn_tris) {
  // fetch features stored in .geogram attributes
  GEO::Attribute<int> attr_corners, attr_se, attr_cce;
  attr_corners.bind_if_is_defined(sf_mesh.vertices.attributes(), "corner");
  attr_se.bind_if_is_defined(sf_mesh.edges.attributes(), "se");
  attr_cce.bind_if_is_defined(sf_mesh.edges.attributes(), "cce");
  if (!attr_corners.is_bound() || !attr_se.is_bound() || !attr_cce.is_bound()) {
    logger().error("ERROR: Mesh feature are not binded!");
    log_and_throw("ERROR");
  }
  for (int v = 0; v < sf_mesh.vertices.nb(); v++) {
    if (attr_corners[v] == 0) continue;
    corners.push_back(v);
  }
  for (int e = 0; e < sf_mesh.edges.nb(); e++) {
    int v1 = sf_mesh.edges.vertex(e, 0);
    int v2 = sf_mesh.edges.vertex(e, 1);
    std::array<int, 2> edge = {{v1, v2}};
    std::sort(edge.begin(), edge.end());
    if (attr_se[e] == 1)
      s_edges.push_back(edge);
    else if (attr_cce[e] == 1)
      cc_edges.push_back(edge);
  }

  // init conn_tris, used by MeshWrapper::get_fs_pair_given_edge()
  for (int i = 0; i < input_faces.size(); i++) {
    const auto& f = input_faces[i];
    for (int j = 0; j < 3; j++) {
      conn_tris[input_faces[i][j]].insert(i);
    }
  }
  // init se_ref_fs_pairs and se_normals
  for (const auto& se : s_edges) {
    std::vector<int> n12_f_ids;
    set_intersection(conn_tris[se[0]], conn_tris[se[1]], n12_f_ids);
    if (n12_f_ids.size() != 2)
      log_and_throw("ERROR: we don't know how to handle if not manifold!!");
    int f1 = n12_f_ids[0];
    int f2 = n12_f_ids[1];
    std::array<int, 2> ref_fs_pair = {{f1, f2}};
    std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
    se_ref_fs_pairs.insert(ref_fs_pair);
    Vector3 n1 = get_mesh_facet_normal(sf_mesh, f1);
    Vector3 n2 = get_mesh_facet_normal(sf_mesh, f2);
    se_normals[se] = {{n1, n2}};
  }

  // sanity check of corners
  // connect to at least 3 sharp edges
  // (not including concave edges)
  std::map<int, std::set<int>> neighbor_v;
  for (const auto& e : s_edges) {
    neighbor_v[e[0]].insert(e[1]);
    neighbor_v[e[1]].insert(e[0]);
  }
  for (const int corner : corners) {
    if (neighbor_v.find(corner) == neighbor_v.end()) {
      logger().error("ERROR: corner {} adjacent to 0 sharp edge", corner);
      log_and_throw("ERROR");
    }
    if (neighbor_v.at(corner).size() <= 2) {
      logger().error("ERROR: corner {} adjacent to <= 2 sharp edge: {}", corner,
                     neighbor_v.at(corner));
      log_and_throw("ERROR");
    }
  }

  logger().debug("#concave_edges = {}", cc_edges.size());
  logger().debug("#sharp_edges = {}", s_edges.size());
  logger().debug("se_normals size: {}", se_normals.size());
  logger().debug("se_ref_fs_pairs size: {}", se_ref_fs_pairs.size());
  logger().debug("#corners  = {}", corners.size());
}

// sf_mesh is matching input_vertices & input_faces
void load_mesh_features(
    const GEO::Mesh& sf_mesh, const std::vector<Vector3>& input_vertices,
    const std::vector<Vector3i>& input_faces,
    std::vector<std::array<int, 2>>& s_edges,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals,
    std::set<std::array<int, 2>>& se_ref_fs_pairs,
    std::vector<std::array<int, 2>>& cc_edges, std::vector<int>& corners,
    std::map<int, std::unordered_set<int>>& conn_tris) {
  s_edges.clear();
  se_normals.clear();
  se_ref_fs_pairs.clear();
  cc_edges.clear();
  corners.clear();
  conn_tris.clear();
  load_mesh_features_helper(sf_mesh, input_vertices, input_faces, s_edges,
                            se_normals, se_ref_fs_pairs, cc_edges, corners,
                            conn_tris);
}

void init_input_faces_normals(const GEO::Mesh& input_mesh,
                              std::vector<Vector3>& input_fnormals) {
  input_fnormals.clear();
  input_fnormals.resize(input_mesh.facets.nb());
  for (size_t f = 0; f < input_mesh.facets.nb(); f++) {
    GEO::vec3 n = GEO::Geom::mesh_facet_normal(input_mesh, f);
    input_fnormals[f] << n[0], n[1], n[2];
    input_fnormals[f].normalize();
  }
}

////////////////////////////////////////////////////////////////////////////////
// Mesh Processsing
// Helper functions
////////////////////////////////////////////////////////////////////////////////

////////////////////////
// For a mesh point:
// density: rho =  1 / (r_sample*LFS)^4
//
// Density will be used later for deciding how many
// random sample points added into input triangle
// (deleted) and LFS_rho and FD_rho will be used for CVT as well!!
void bind_density_to_input_mesh_vertices(
    const GEO::NearestNeighborSearch_var& lfs_kd_tree,
    const std::vector<double>& lfs_rho_values, const double& r_sample,
    GEO::Mesh& sf_mesh, bool is_debug = false) {
  if (is_debug) logger().debug("start binding mesh vertex weight/density ...");
  // sanity checks
  if (lfs_kd_tree == nullptr || lfs_rho_values.empty()) {
    log_and_throw("ERROR: LFS must be initialized first!!!");
  }
  if (lfs_rho_values.size() != lfs_kd_tree->nb_points()) {
    logger().debug("lfs_kd_tree nb_points: {}", lfs_kd_tree->nb_points());
    logger().debug("rho_max_values size: {}", lfs_rho_values.size());
    log_and_throw("ERROR: rho_max_values.size() != lfs_kd_tree.nb_points()");
  }

  // bind weight/density info to mesh vertices
  // density: rho = 1/ (r*LFS)^4
  GEO::Attribute<double> densities;
  densities.bind(sf_mesh.facets.attributes(), "density");

  GEO::index_t K_nb_neigh = 1;  // k nearest neighbors
  std::vector<GEO::index_t> k_neighbors(K_nb_neigh);
  std::vector<double> k_neighbors_sq_dist(K_nb_neigh);

  // LFS
  // we store nearest max lfs_rho for each point, no matter feature or not
  for (int f = 0; f < sf_mesh.facets.nb(); f++) {
    densities[f] = -1;  // init value
    for (int lv = 0; lv < sf_mesh.facets.nb_vertices(f); lv++) {
      const int vidx = sf_mesh.facets.vertex(f, lv);
      const GEO::vec3& p_geo = sf_mesh.vertices.point(vidx);
      std::vector<double> p = {{p_geo[0], p_geo[1], p_geo[2]}};

      lfs_kd_tree->get_nearest_neighbors(
          K_nb_neigh, p.data(), k_neighbors.data(), k_neighbors_sq_dist.data());
      double max_lfs_rho = DBL_MIN;
      for (const GEO::index_t kn_id : k_neighbors) {
        max_lfs_rho = std::max(max_lfs_rho, lfs_rho_values[kn_id]);
      }
      densities[f] = std::max(max_lfs_rho, densities[f]);
    }  // for lv of facet
  }    // for shape3D.sf_mesh.facets
}

////////////////////////////////////////////////////////////////////////////////
// Mesh Processsing
// Main functions
////////////////////////////////////////////////////////////////////////////////
void mesh_remesh_split_sharp_edges(ThreeDimensionalShape& shape3D,
                                   bool is_debug) {
  if (is_debug) logger().debug("calling mesh_remesh_split_sharp_edges ...");
  // bind density to each LFS mesh vertex
  bind_density_to_input_mesh_vertices(
      shape3D.lfs_kd_tree, shape3D.rho_max_values, shape3D.params.r_sample,
      shape3D.sf_mesh, is_debug);
  ///////////////////////////////////////////////////////////////////////////////////////
  const GEO::Mesh& mesh = shape3D.sf_mesh;
  const MeshWrapper& sf_mesh_wrapper = shape3D.sf_mesh_wrapper;
  const AABBWrapper& aabb_wrapper = shape3D.aabb_wrapper;
  const Parameters& params = shape3D.params;

  std::vector<double>& seed_points = shape3D.seed_points;
  std::vector<Vector3>& sf_seeds = shape3D.sf_seeds;
  std::vector<Vector3>& seed_normals_v2 = shape3D.seed_normals_v2;
  std::set<int>& feature_points = shape3D.feature_points;
  std::set<int>& point_is_locked = shape3D.point_is_locked;
  seed_points.clear();
  sf_seeds.clear();
  feature_points.clear();
  point_is_locked.clear();
  // TODO: estimated density might < original mesh density
  size_t nb_points = params.init_nb_samples;
  if (is_debug) logger().info("Ideal #samples: {}", nb_points);

  // LFS
  const GEO::NearestNeighborSearch_var& lfs_kd_tree = shape3D.lfs_kd_tree;
  const std::vector<double>& lfs_min_values = shape3D.lfs_min_values;
  const double s_edge_ideal_eps_rel = shape3D.s_edge_ideal_eps_rel;

  const double& ideal_length = params.init_edge_length;
  const int dim = 3;
  if (nb_points < mesh.vertices.nb()) {
    nb_points = mesh.vertices.nb();
    logger().info("Current model meets our desired density, keep #point {} ",
                  nb_points);
  }

  // currently all indices mathcing sf_mesh
  // but will be updated to matching sf_seeds/seed_points
  std::vector<int>& corners = shape3D.corners;
  std::set<int> corners_set;  // for fast indexing
  for (const auto& c : corners) corners_set.insert(c);
  std::vector<std::array<int, 2>>& s_edges = shape3D.s_edges;
  std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals =
      shape3D.se_normals;
  std::vector<std::array<int, 2>>& cc_edges = shape3D.cc_edges;

  // TODO: replace corners, sharp edges, concave edges
  std::vector<SpecialCorner>& sp_corners = shape3D.sp_corners;  // corners
  std::vector<SpecialEdge>& sp_edges =
      shape3D.sp_edges;  // sharp edges and concave edges
  std::map<std::array<int, 2>, int>& map_to_sp_edges = shape3D.map_to_sp_edges;
  std::vector<TangentConcaveLine>& tan_cc_lines = shape3D.tan_cc_lines;

  // temporary variables
  std::map<int, int> seeds_map;  // indices from sf_mesh to sf_seeds/seed_points
  std::vector<int> corners_new;
  std::map<int, Vector3> corners_normals_new;
  std::vector<std::array<int, 2>> s_edges_new;
  std::map<std::array<int, 2>, std::array<Vector3, 2>> se_normals_new;
  std::vector<std::array<int, 2>> cc_edges_new;

  if (is_debug) logger().info("start locking sharp edges ...");
  load_sample_and_lock_sharp_edges(
      sf_mesh_wrapper, ideal_length, s_edge_ideal_eps_rel, lfs_kd_tree,
      lfs_min_values, params.r_sample, corners_set, seed_points, seeds_map,
      point_is_locked, s_edges, se_normals, sp_edges, map_to_sp_edges,
      is_debug);
  shape3D.reload_s_edges_adjs();

  if (is_debug) logger().info("start loading concave edges (no lock) ...");
  // will init tan_cc_lines
  load_concave_edges(sf_mesh_wrapper, seed_points, seeds_map, cc_edges,
                     sp_edges, tan_cc_lines, map_to_sp_edges, is_debug);
  shape3D.reload_ref_fs_pairs_not_cross();

  if (is_debug) logger().info("start locking corners ...");
  load_and_lock_corners(sf_mesh_wrapper, seed_points, seeds_map,
                        point_is_locked, corners, is_debug);
  // NOTE: corners are mapping to seed_points now
  load_corner_special_edges(sp_edges, corners, sp_corners, is_debug);
  if (is_debug)
    logger().debug("#seeds {}/{} added", seed_points.size() / 3, nb_points);

  if (nb_points < seed_points.size() / 3) {
    int new_nb_points = nb_points + seed_points.size() / 3;
    if (is_debug)
      logger().info(
          "we need more samples to keep feature, expand #seeds: {} -> {}",
          nb_points, new_nb_points);
    nb_points = new_nb_points;
  }

  //////////
  // TODO: deprecate this, we update seed normals after CVT projection
  // Update seed_normals
  // untill now, we only have feature seeds, therefore, insert 0 to
  // seed_normals feature seeds normals are stored in se_normals
  seed_normals_v2.clear();
  seed_normals_v2.resize(seed_points.size() / 3, Vector3(0, 0, 0));

  // reload feature points from corners & sharp edges
  // NOTE: feature_points != point_is_locked
  shape3D.reload_feature_points();

  //////////
  // Insert random samples
  // = init_nb_samples - locked sample size
  int nb_new_points = nb_points - point_is_locked.size();
  GEO::Attribute<double> _;
  // random surface seeds
  std::vector<double> points_new;
  std::vector<double> points_new_normals;

  if (is_debug) logger().debug("going to add random points manually ...");

  if (params.downsample_percentage != -1) {
    if (is_debug) {
      logger().debug(
          "given downsample_percentage: {}, please make sure samples are DENSE "
          "enough!",
          params.downsample_percentage);
    }
    mesh_generate_adaptive_samples_on_surface_with_normals_given_downsample<
        dim>(mesh, points_new, points_new_normals,
             params.downsample_percentage);
  } else {
    if (is_debug) logger().debug("going to add random points using LFS ...");
    // here we don't care about nb_new_points,
    // allow the density to handle the number of seeds
    GEO::Attribute<double> facet_density_attr(mesh.facets.attributes(),
                                              "density");
    mesh_generate_adaptive_samples_on_surface_with_normals_given_sizing_field<
        dim>(mesh, points_new, points_new_normals, facet_density_attr);
  }

  nb_new_points = points_new.size() / 3;
  if (is_debug) {
    logger().info("added {} random samples", nb_new_points);
    logger().debug("points_new.size(): {}", points_new.size() / 3);
    logger().debug("points_new_normals.size(): {}",
                   points_new_normals.size() / 3);
  }

  // save points to seed_points and seed_normals_v2
  for (int i = 0; i < points_new.size() / 3; i++) {
    // store seeds and seeds normals
    for (int j = 0; j < 3; j++) {
      seed_points.push_back(points_new[i * 3 + j]);
    }
    seed_normals_v2.push_back(Vector3(points_new_normals[i * 3],
                                      points_new_normals[i * 3 + 1],
                                      points_new_normals[i * 3 + 2]));
  }
  nb_points = seed_points.size() / dim;
  if (is_debug)
    logger().debug("#all seed_points size {}, seed_normals_v2: {}", nb_points,
                   seed_normals_v2.size());

  // sanity checks
  if (nb_points != seed_normals_v2.size()) {
    logger().error("nb_points: {} but seed_normals_v2 {}", nb_points,
                   seed_normals_v2.size());
    log_and_throw("ERROR");
  }
  if (seed_points.size() % dim != 0) {
    log_and_throw("ERROR: seed_points.size() must be mutiplication of 3");
  }

  if (is_debug)
    logger().debug("#seeds {}/{} added", seed_points.size() / dim, nb_points);
  if (nb_points * 3 != seed_points.size()) nb_points = seed_points.size() / dim;

  // copy shape3D.seed_points to shape3D.init_seed_points
  shape3D.init_seed_points.clear();
  std::copy(seed_points.begin(), seed_points.end(),
            std::back_inserter(shape3D.init_seed_points));

  // TODO: replace seed_points by sf_seeds
  // copy to sf_seeds
  for (int i = 0; i < seed_points.size() / 3; i++) {
    sf_seeds.push_back(Vector3(seed_points[i * 3], seed_points[i * 3 + 1],
                               seed_points[i * 3 + 2]));
  }
  if (is_debug) logger().debug("done mesh_remesh_split_sharp_edges");
}

}  // namespace matfp