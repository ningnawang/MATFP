// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/delaunay/delaunay.h>
#include <geogram/points/kd_tree.h>
#include <geogram/points/nn_search.h>

#include <Eigen/Dense>

#include "matfp/AABBWrapper.h"
#include "matfp/Args.h"
#include "matfp/Logger.h"
#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/Parameters.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/Types/RegularTriangulation.h"

namespace matfp {

enum EdgeType {
  UE = -1,  // unkown edge
  SE = 1,   // convex sharp edge
  CCE = 2   // concave edge
};

// can be sharp edge or concave edge
class SpecialEdge {
 public:
  SpecialEdge(const int _id, const EdgeType& _t,
              const std::array<int, 2>& _vs_ids);
  ~SpecialEdge(){};

  void print_info() const;

  int id;
  EdgeType type;
  std::array<int, 2> vs_seed_ids;  // mapping to seed_points
  std::array<Vector3, 2> vs_pos;
  std::array<int, 2> adj_ref_fs_pair;         // mapping to GEO::Mesh sf_mesh
  std::array<Vector3, 2> adj_ref_tan_points;  // centroid of adj_ref_fs_pair
  std::array<Vector3, 2> adj_ref_normals;     // order maps adj_ref_fs_pair

  // TODO: remove this?
  // if this is a concave edge
  // we need to know which TangentConcaveLine it is mapping to
  int cc_line_id = -1;
  // mapping to MVertex::tag
  // could be {-1,-1} for concave edges
  std::array<int, 2> mat_se_tags = {{-1, -1}};

  // for corner preservations
  // only sharp edge has id != -1
  // concave line have id = -1
  int sheet_node_id = -1;
};

class SpecialCorner {
 public:
  SpecialCorner(const int _id, const int _s_id, const Vector3& _pos);
  ~SpecialCorner(){};

  void print_info() const;
  void push_to_adj_edges(const int sp_eid);

  int id;
  int seed_id;  // mapping to seed_points
  Vector3 pos;
  std::set<int> adj_edges;            // mapping to SpecialEdge
  std::vector<int> adj_edges_sorted;  // mapping to SpecialEdge
  int mat_tag = -1;             // each corner must be a zero-radius sphere
  double min_corner_len = DBL_MAX;  // define a small area around corner, updated in
                                // func load_corners_min_length()
};

// wrapper of GEO::Mesh sf_mesh
class MeshWrapper {
 public:
  Vector3 get_v_pos(const int vid) const;
  void get_fs_pair_given_edge(const std::array<int, 2>& e,
                              std::vector<int>& fs_pair) const;
  Vector3 get_f_normal(const int fid) const;
  Vector3 get_f_centroid(const int fid) const;

  std::vector<Vector3> input_vertices;
  std::vector<Vector3i> input_faces;
  std::vector<Vector3> input_fnormals;  // face normals
  std::vector<int> input_tags;
  std::map<int, std::unordered_set<int>> conn_tris;  // map from v -> fs

  Eigen::MatrixXd VI;
  Eigen::MatrixXi FI;
};

class ThreeDimensionalShape {
 public:
  Parameters params;
  std::string mesh_name;

  NonManifoldMesh mat_refined;

  //////////////////////////////////
  // Input Mesh (after subdivision)
  //////////////////////////////////
  GEO::Mesh sf_mesh;  // mesh after subdivision
  // new representation of sf_mesh
  MeshWrapper sf_mesh_wrapper;

  /////////////////////////////////////////////
  // Sharp edges (may differ from input mesh)
  // mapping to seed_points
  /////////////////////////////////////////////
  std::vector<SpecialEdge> sp_edges;  // TODO: replace s_edges and cc_edges
  std::map<std::array<int, 2>, int>
      map_to_sp_edges;                    // from s_edges/cc_edges to sp_edges
  std::vector<SpecialCorner> sp_corners;  // TODO: replace corners

  std::vector<std::array<int, 2>> s_edges;  // sharp edges (se)
  std::map<std::array<int, 2>, std::array<Vector3, 2>>
      se_normals;  // matches s_edges, not GEO::Mesh sf_mesh
  std::set<std::array<int, 2>>
      se_ref_fs_pairs;  // matches GEO::Mesh sf_mesh, not s_edges,
                        // used by iterate_sphere()
  std::set<std::array<int, 2>>
      ref_fs_pairs_not_cross;  // se_ref_fs_pairs + tan_cc_lines.adj_ref_fs_pair
  std::vector<std::array<int, 2>>
      cc_edges;  // concave edges, if cc, then not se, may be updated
  // id mapping to cc_edges after subdivision
  // but cc_edges then will be updated to seed_points
  std::vector<int> corners;  // store vertex indices
  std::map<int, std::set<int>>
      neighbors;  // [no use] 1-ring neighbors of vertices
  const double s_edge_ideal_eps_rel = 1;  // (in % of the sizing field)

  GEO::Mesh c_mesh;  // corner mesh for special sampling around corners, not
                     // controlled by sizing field
  std::vector<std::array<int, 3>> c_faces;  // use for update c_mesh
  std::vector<Vector3> c_mid_vertices;  // this store middle vertices created
                                        // between two sharp edge along corners
  std::vector<Vector3> c_mid_normals;   // normals for c_mid_vertices

  // Local Feature Size
  GEO::NearestNeighborSearch_var lfs_kd_tree;
  std::vector<double> lfs_seeds;       // non-feature seeds points (= dt_vs =
                                       // seed_points - feature_points)
  std::vector<double> lfs_min_values;  // size = lfs_seeds/3
                                       // TODO: rename it to lfs_rho_values
  std::vector<double>
      rho_max_values;  // rho = 1. / std::pow(r_sample * minLFS, 4)
  double max_lfs_value;

  std::vector<double> init_seed_points;  // all seeds before CVT
  std::vector<double>
      seed_points;  // all seeds [deprecate? use sf_seeds instead?]
  std::vector<Vector3> sf_seeds;  // replacing seed_points
  std::set<int> seed_is_deleted;  // delete seed after CVT if inside c_mesh
  // store each seed point's normal regarding input facets
  //
  // feature seeds -> (0, 0, 0), please use se_normals or corners_normals
  // non-feature seeds -> (x, y, z)
  std::vector<Vector3> seed_normals_v2;  // replace seed_normals
  std::set<int> feature_points;   // a set of indices indicate the locked seeds,
                                  // regarding this->seed_points
  std::set<int> point_is_locked;  // includes feature_points and other locked
                                  // points (vertices created by c_mesh)

  // Delaunay Triangulation
  DelaunayTriangulation dt;

  ///////////
  // for Regular Triangulation / Weighted Delaunay Triangulation
  // using power distance
  // TODO:
  // for replacing all below nonsense
  std::vector<MVertex> all_medial_spheres;

  // init the same as cc_edges, but store more
  // all concave lines from SpecialEdge of type EdgeType::CCE
  std::vector<TangentConcaveLine> tan_cc_lines;

  // map from rt tag to all_medial_spheres
  // rt tag only store valid (!is_outside && !is_deleted) medial spheres
  // to reduce the complexity of RPD and MAT creation later on
  // this can also be called map_rt_tag_to_all, which one you like?
  std::vector<int> valid_medial_spheres;
  std::map<int, int> all_to_valid_medial_spheres;  // deprecating
  std::set<int> sphere_for_debug;  // tag here is mapping all_medial_spheres,
                                   // not valid_medial_spheres

  // sharp edge kd-tree
  // save se_spheres in another way
  // map of a sharp edge seed and its two connected sharp edge seeds
  std::map<int, std::unordered_set<int>>
      s_edges_adjs;  // deprecating, using se_spheres_adjs instead
  std::set<std::array<int, 2>>
      se_spheres;  // store sharp edge pairs in MVertex::tag
  // std::map<int, std::unordered_set<int>> se_spheres_adjs;  // in MVertex::tag
  GEO::NearestNeighborSearch_var se_kd_tree;
  std::vector<double> se_kd_points;  // this storing is necessary, otherwise
                                     // se_kd_tree dunno what to search
  std::map<int, int>
      se_kd_tree_idx_to_se_tag;  // idx from se_kd_points to MVertex::tag
  void reload_s_edges_adjs();

  RegularTriangulationNN rt;

  std::vector<double> winding_num;  // non-feature mat pts only
  std::vector<bool> is_outside;     // non-feature mat pts only
  std::vector<double>
      all_weighted_mat_pts_cgal;  // exclude outside non-feature mat pts, dim=4
  std::vector<bool> is_locking;   // index over all_weighted_mat_pts_cgal

  // for RPD (power RVD)
  bool is_volumetric = false;
  GEO::Mesh power_rvd;    // geogram
  GEO::Mesh rpd_refined;  // RPD with refinement
  GEO::Mesh
      rpd_merged;  // RPD with only simple mesh representation of powercell
  Vertex_handle_rt current_seed_handle;  // for debug only
  int current_facet;                     // for debug use only
  // RPD surface points and their adjacent points
  // rpd_seed_adj is shared by both rpd_no_refinement and rpd_refined
  // indices are matching valid_medial_spheres, not all_medial_spheres
  std::map<GEO::index_t, std::set<GEO::index_t>> rpd_seed_adj;
  std::map<GEO::index_t, std::set<GEO::index_t>>
      rpd_vs_bisectors;  // each vertex of RPD, store its bisectors
  // used to update RPD mesh, merge vertices that are close
  // and merge vertices in degenerated faces
  // std::map<GEO::index_t, GEO::index_t> rpd_v_old2new;
  std::map<GEO::index_t, std::set<std::array<GEO::index_t, 2>>>
      rpd_cell_boundary_segments;

  AABBWrapper aabb_wrapper;
  Args args;  // deprecating

 public:
  void reload_ref_fs_pairs_not_cross();

 public:
  //// initialization
  inline Scalar get_sf_diag() const { return GEO::bbox_diagonal(sf_mesh); }

  inline void reload_feature_points() {
    // reload feature points from corners and sharp edges
    feature_points.clear();
    for (const auto& c : corners) {
      feature_points.insert(c);
    }
    for (auto const& se_pair : s_edges) {
      feature_points.insert(se_pair[1]);
      feature_points.insert(se_pair[0]);
    }
  }

  // DEBUG ONLY
  // loop valid_medial_spheres to find valid idx
  // if not found, we just return -1
  inline int find_valid_medial_sphere_given_all_idx(const int all_idx) {
    int valid_idx = -1;
    for (int i = 0; i < valid_medial_spheres.size(); i++) {
      if (valid_medial_spheres[i] == all_idx) {
        valid_idx = i;
        break;
      }
    }
    return valid_idx;
  }

  inline void set_bb_corners() {
    Vector3 min, max;
    min = sf_mesh_wrapper.input_vertices.front();
    max = sf_mesh_wrapper.input_vertices.front();

    for (size_t j = 0; j < sf_mesh_wrapper.input_vertices.size(); j++) {
      for (int i = 0; i < 3; i++) {
        min(i) = std::min(min(i), sf_mesh_wrapper.input_vertices[j](i));
        max(i) = std::max(max(i), sf_mesh_wrapper.input_vertices[j](i));
      }
    }
    // const Scalar dis = std::max(params.init_edge_length, params.eps_input *
    // 2); for (int j = 0; j < 3; j++) { 	min[j] -= dis;xmin
    params.bbox_min = min;
    params.bbox_max = max;

    logger().debug("min = {} {} {}", min[0], min[1], min[2]);
    logger().debug("max = {} {} {}", max[0], max[1], max[2]);
  }

  inline void set_bb_corners(GEO::Mesh& mesh) {
    GEO::vec3& p0 = mesh.vertices.point(0);
    Vector3 min(p0[0], p0[1], p0[2]);
    Vector3 max(p0[0], p0[1], p0[2]);

    for (int i = 0; i < (int)mesh.vertices.nb(); ++i) {
      GEO::vec3& p = mesh.vertices.point(i);
      for (int j = 0; j < 3; j++) {
        min(j) = std::min(min(j), p[j]);
        max(j) = std::max(max(j), p[j]);
      }
    }
    params.bbox_min = min;
    params.bbox_max = max;

    logger().debug("min = {} {} {}", min[0], min[1], min[2]);
    logger().debug("max = {} {} {}", max[0], max[1], max[2]);
  }

  inline void set_bb_corners(NonManifoldMesh& mat) {
    mat.computebb();
    params.bbox_min = Vector3(mat.m_min[0], mat.m_min[1], mat.m_min[2]);
    params.bbox_max = Vector3(mat.m_max[0], mat.m_max[1], mat.m_max[2]);
  }

 public:
  ThreeDimensionalShape(){};
  ~ThreeDimensionalShape() {
    se_kd_tree.reset();
    lfs_kd_tree.reset();
    rt.clear();
  };
};

}  // namespace matfp