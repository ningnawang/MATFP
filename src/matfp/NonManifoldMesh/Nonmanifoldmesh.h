// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <float.h>
#include <geogram/points/kd_tree.h>
#include <geogram/points/nn_search.h>
#include <matfp/Logger.h>
#include <matfp/Types/CommonTypes.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <unordered_set>
#include <vector>

#include "Primitives.h"
#include "matfp/LinearAlgebra/Wm4Matrix.h"
#include "matfp/LinearAlgebra/Wm4Vector.h"

namespace matfp {

using namespace std;

class NonManifoldMesh_Vertex {
 public:
  std::set<int> edges_;  // edges that connect with neighbor mat vertices
  std::set<int> faces_;  // mat faces that contain this vertex

  bool HasEdge(int eid) { return (edges_.find(eid) != edges_.end()); }
  bool HasFace(int fid) { return (faces_.find(fid) != faces_.end()); }

 public:
  // index of inside mat vertices,
  // same as RT::RVI::tag, given when calculating MAT
  // // tag is matching index of all_weighted_mat_pts_cgal
  // tag is matching index of valid_medial_spheres
  //
  // ninwang: 2022-02-21
  // may not matching valid_idx
  int tag;
  int all_idx;  // matching MVertex::tag
  SphereType type = SphereType::T_UNK;

  double winding_number = -1;
  bool is_outside = false;
  bool is_deleted = false;

  // TODO: start using these
  Vector3 pos;
  Vector2 pos_2d;
  double radius = 0;  // 0 is when dim=3, no weight info

 public:
  bool operator<(NonManifoldMesh_Vertex const& b) { return tag < b.tag; }
  static bool comp_by_id(const NonManifoldMesh_Vertex* a,
                         const NonManifoldMesh_Vertex* b) {
    return a->tag < b->tag;
  };
  void print_mat_vertex();
};

class NonManifoldMesh_Edge {
 public:
  std::pair<int, int> vertices_;  // vertex list
  std::set<int> faces_;           // triangle list
  bool is_deleted = false;

  std::array<int, 2> all_vs_tags;
  int rei_id;  // for debug
  std::vector<Vector3> dual_polygon;
  std::vector<bool> is_dual_poly_vs_outside;
  // 1. any dual vs is inside
  // 2. dual poly has intersect
  bool is_dual_poly_intersect = false;
  std::vector<int> poly_bbox_isect_triangles;

 public:
  void print_mat_edge() const;

  bool HasVertex(int vid) {
    return ((vertices_.first == vid) || (vertices_.second == vid));
  }
  bool HasFace(int fid) { return (faces_.find(fid) != faces_.end()); }

 public:
  NonManifoldMesh_Edge() {}

 public:
  int tag;
  double length;

  Cone cone;
  bool is_valid_cone;  // true when type != 1 (height != 0.0)
};

class NonManifoldMesh_Face {
 public:
  // std::vector<int> vertices_sorted_; // when order matters

  // set has its own order
  std::set<int> vertices_;  // vertex list
  std::set<int> edges_;     // edge list
  std::set<int> tets_;      // neighboring tetrahedrons
  std::set<int> tets_old_;  // tets before thinning
  bool is_deleted = false;

  bool HasVertex(int vid) { return (vertices_.find(vid) != vertices_.end()); }
  bool HasEdge(int eid) { return (edges_.find(eid) != edges_.end()); }

  // these are called after NonManifoldMesh is created
  // updated by calling ComputeFacesCentroid()
  Vector3 centroid;

  // updated by calling ComputeFacesNormal()
  // note: this is not a 'real' normal, just normal of given orientation for
  // rendering
  Vector3 normal;

  std::array<int, 3> all_vs_tags = {{-1, -1, -1}};
  int rfi_id;  // for debug
  // the dual segment of RT face
  // assigned by function generate_RT_dual_info()
  std::array<Vector3, 2> dual_segment;
  std::array<bool, 2> is_dual_seg_vs_outside;
  bool is_dual_seg_intersect = false;
  // intersection between dual_segment and surface
  std::array<Vector3, 2> dual_intersections;
  double dist_dual_intersections = -1;

  // importance metric
  bool is_locked = false;
  double importance = 2.;
  int delete_for_tet = -1;

 public:
  void print_mat_face() const;

 public:
  enum Type {
    INVALID = 1,      // -- fail to generate 2 simple triangles
    TETRAHEDRON = 2,  // -- share an edge
    PRISM = 3         // -- share a point or nothing
  };

 public:
  int tag;
  Type type;

 public:
  //////////////////////////////////////////////////////////////////////
  // There are 3 vertex pairs, each pairs share the same sphere
  // (must be on the opposite side of a same sphere)
  // the pairs are:
  // st0.v[0] - st1.v[0]
  // st0.v[1] - st1.v[2]
  // st0.v[2] - st1.v[1]
  //
  // Also 3 vertices in each SimpleTriangle are counter-clockwise
  //
  // Detail see Primitives::get_triangles_from_three_spheres()
  //////////////////////////////////////////////////////////////////////
  SimpleTriangle st[2];
  bool is_valid_st;
  double st_hausdorff_dist[2][3];
};

// These tets need to be removed
class NonManifoldMesh_Tet {
 public:
  int tag;
  // set has its own order
  std::set<int> vertices_;    // vertex list
  std::set<int> edges_;       // edge list
  std::set<int> faces_;       // facet list
  std::set<int> neigh_tets_;  // neighboring tetrahedrons
  bool is_deleted = false;

  Vector3 dual_vertex;
  bool is_dual_vertex_outside = false;
  double dual_sq_radius;
  std::array<int, 4> all_vs_tags;

  bool HasVertex(int vid) { return (vertices_.find(vid) != vertices_.end()); }
  bool HasEdge(int eid) { return (edges_.find(eid) != edges_.end()); }
  bool HasFace(int fid) { return (faces_.find(fid) != faces_.end()); }

 public:
  void print_mat_tet() const;
};

typedef std::pair<bool, NonManifoldMesh_Vertex*> Bool_VertexPointer;
typedef std::pair<bool, NonManifoldMesh_Edge*> Bool_EdgePointer;
typedef std::pair<bool, NonManifoldMesh_Face*> Bool_FacePointer;

class NonManifoldMesh {
 public:
  NonManifoldMesh() {
    seconds = 0.;
    numVertices = numEdges = numFaces = numTets = 0;
    numVertices_active = numEdges_active = numFaces_active = numTets_active = 0;
    vs_adj_pair_to_skip.clear();
  }

 public:
  int numVertices;         // including both feature and non-feature points
  int numVertices_active;  // no use
  int numEdges;
  int numEdges_active;
  int numFaces;
  int numFaces_active;
  int numTets;
  int numTets_active;

 public:
  double max_radius = DBL_MIN;

  std::map<int, int> eid2rei;  // for debug, NonManifoldMesh_Edge::id -> REI::id

 public:
  // store mat edges that should not exist
  // this should mainly happend around concave lines
  // updated by function add_vs_pair_to_skip()
  std::set<std::array<int, 2>> vs_adj_pair_to_skip;

  // check the invalid connection (edge) of MAT
  std::set<std::array<int, 2>> invalid_mat_edges;

  //////////// Unthin !!!
  // check the thinness of MAT
  // indices matches valid_medial_spheres
  std::map<int, std::unordered_set<int>> mat_v_neighbors;
  std::set<std::vector<int>> tets_in_mat;
  void check_and_store_unthin_tets_in_mat(
      bool is_force_update_neighs = false);  // update tets_in_mat

  // problematic faces that makes euler increase
  std::set<int> unthin_faces;
  std::vector<int> unthin_trace;
  std::queue<int> next_fs_to_check;
  std::set<int> vs_visited;  // mat vertices
  std::set<int> es_visited;  // mat edges
  std::set<int> fs_visited;  // mat faces
  int num_undeleted_fs = 0;  // < fs_visited.size()
  int prev_euler = 1;
  void trace_unthin_mat_faces();

 public:
  // for external medial features
  std::set<std::array<int, 2>> mat_extf_edges;

  // for internal medial features
  std::set<std::array<int, 2>> mat_intf_edges;
  std::map<int, std::unordered_set<int>> mat_intf_adjs;  // valid_idx
  bool is_feature_edge_exists(const int v1, const int v2) const;
  void insert_new_feature_edge(const int v1, const int v2);

  // TODO: deprecating
  // indices mapping to valid_idx
  GEO::NearestNeighborSearch_var intf_kd_tree;
  std::vector<double> intf_kd_points;  // this storing is necessary, otherwise
                                       // intf_kd_tree dunno what to search
  std::map<int, int> intf_kd_tree_idx_to_intf_adjs;  // idx from intf_kd_points
                                                     // to mat_intf_edges key
  void reload_intf_kd_tree();

 public:
  // in case people import .ma directly
  double m_min[3];
  double m_max[3];
  double bbox_diag_l;

  double seconds;
  std::string mat_name;

 public:
  std::vector<Bool_VertexPointer> vertices;
  std::vector<Bool_EdgePointer> edges;
  std::vector<Bool_FacePointer> faces;
  std::vector<NonManifoldMesh_Tet> tets;  // the one to remove!!!

 public:
  void print_info();
  bool is_empty();
  int get_sphere_all_idx(const int& valid_idx);

  int create_vertex(const double x, const double y, const double z,
                    const double radius,
                    int tag = -1 /*if not given then auto-increase*/,
                    const bool is_outside = false,
                    const SphereType type = SphereType::T_UNK,
                    const int all_idx = -1);
  int create_edge(const int vid0, const int vid1,
                  int tag = -1 /*if not given then auto-increase*/);
  int create_face(const std::set<int>& vvid,
                  int tag = -1 /*if not given then auto-increase*/);
  int create_tet(const std::vector<int> fids,
                 int tag = -1 /*if not given then auto-increase*/,
                 bool is_auto_create_face = false /*create face if not exist*/);

  void clear();
  int clear_all_tets();

  bool delete_edge(const int eid);
  bool delete_face(const int fid);
  bool delete_tet(const int tet_id);

  void get_adj_vertices(int vid, std::unordered_set<int>& adj_vs) const;

  // use vs_adj_pair_to_skip to skip creating
  // mat edges or faces
  // NOTE: please load vs_adj_pair_to_skip before using
  void add_vs_pair_to_skip(int v1, int v2);
  bool is_skip_mat_edge_creation(int v1, int v2);
  bool is_skip_mat_face_creation(const std::set<int>& vvid);

  // given a point from input mesh,
  // return the nearest distance from mat,
  // by looping all mat vertices & edges & facets
  // this is mostly used for caculating hausdorff_distance
  double get_nearest_distance(Vector3& input_point);

 public:
  void AdjustStorage();

 public:
  bool ValidVertex(int vid);
  bool Edge(int vid0, int vid1, int& eid);
  bool Face(const std::set<int>& vset, int& fid);
  void UpdateCentroid(int fid);
  void ComputeFacesCentroid();
  void UpdateNormal(int fid);
  void ComputeFacesNormal();
  void GetNeighborVertices(int vid, std::set<int>& neighborvertices);
  void GetLinkedEdges(int eid, std::set<int>& neighboredges);
  void GetAdjacentFaces(int fid, std::set<int>& neighborfaces);
  bool Contractible(int vid_src, int vid_tgt);
  bool MergeVertices(int vid_src1, int vid_src2, int& vid_tgt);

 public:
  void Export(std::string fname);

  void DeleteFace(int fid);
  void DeleteEdge(int eid);
  void DeleteVertex(int vid);

  void InsertVertex(NonManifoldMesh_Vertex* vertex, int& vid);
  void InsertEdge(int vid0, int vid1, int& eid);
  void InsertFace(std::set<int> vset);

  int VertexIncidentEdgeCount(int vid);
  int VertexIncidentFaceCount(int vid);
  int EdgeIncidentFaceCount(int eid);

 public:
  void computebb();

 public:
  void compute_edges_cone();
  void compute_faces_simple_triangles();

 private:
  void compute_edge_cone(int eid);
  void compute_face_simple_triangles(int fid);

};  // NonManifoldMesh

class MatIO {
 private:
  static void get_mat_clean(const NonManifoldMesh& mat,
                            std::vector<Vector4>& vertices,
                            std::vector<std::array<int, 2>>& edges,
                            std::vector<std::array<int, 3>>& faces);
  static void export_ma_given(const std::string& maname,
                              const std::vector<Vector4>& mat_vertices,
                              const std::vector<std::array<int, 2>>& mat_edges,
                              const std::vector<std::array<int, 3>>& mat_faces,
                              bool is_use_given_name = false);

 public:
  // load the user defined mat
  static void load_nmm(const std::string& path, NonManifoldMesh& mat);

  static void export_ma_clean(const std::string& maname,
                              const NonManifoldMesh& mat);

  static void export_nmm(const std::string& maname, const NonManifoldMesh& mat);
  // save as .r file (for ET)
  static void export_nmm_vs(const std::string& maname,
                            const NonManifoldMesh& mat);

  static void write_nmm_ply(const std::string& maname,
                            const NonManifoldMesh& mat);

  static void export_nmm_tets(const std::string& maname,
                              const NonManifoldMesh& mat);
};  // MatIO

}  // namespace matfp
