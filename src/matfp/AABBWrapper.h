// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_geometry.h>

#include <memory>

#include "Types/CommonTypes.h"

namespace matfp {

class AABBWrapper {
 public:
  // GEO::Mesh sf_mesh;  // a copy of real sf_mesh before subdivision
  GEO::Mesh s_mesh;  // [no use] sharp feature
  bool is_s_mesh_exist = false;
  GEO::Mesh c_mesh;   // corner sphere restricted on sf_mesh, this is a copy
  GEO::Mesh cc_mesh;  // for concave lines
  bool is_cc_mesh_exist = false;

 private:
  std::shared_ptr<GEO::MeshFacetsAABB> s_tree;  // sharp edge tree, no use
  std::shared_ptr<GEO::MeshFacetsAABB> sf_tree;
  std::shared_ptr<GEO::MeshFacetsAABB> c_tree;  // no use
  std::shared_ptr<GEO::MeshFacetsAABB> cc_tree;

 public:
  AABBWrapper() {}

  // is_reorder == true => will reorder sf_mesh
  void init_sf_mesh_and_tree(GEO::Mesh &_sf_mesh, bool is_reorder = true) {
    sf_tree = std::make_shared<GEO::MeshFacetsAABB>(_sf_mesh, is_reorder);
  }

  void init_feature_meshes_and_trees(
      const std::vector<Vector3> &input_vertices,
      const std::vector<std::array<int, 2>> &s_edges,
      const std::vector<std::array<int, 2>> &cc_edges);

 private:
  bool init_mesh_from_edges(const std::vector<Vector3> &input_vertices,
                            const std::vector<std::array<int, 2>> &edges,
                            GEO::Mesh &mesh);

 public:
  inline bool is_corner_mesh_valid() const { return c_tree != nullptr; }

  inline Scalar get_sq_dist_to_c_mesh(const Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    c_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    // logger().debug("point ({},{},{}) has dist {} to c_mesh", geo_p[0],
    // geo_p[1], geo_p[2], sq_dist);
    return sq_dist;
  }

  inline bool is_close_to_c_mesh(const Vector3 &p, const Scalar eps_2,
                                 double &sq_dist) const {
    sq_dist = get_sq_dist_to_c_mesh(p);
    if (Scalar(sq_dist) < eps_2) return true;
    return false;
  }

 public:
  inline bool sf_segment_intersection(const Vector3 &p1,
                                      const Vector3 &p2) const {
    GEO::vec3 q1(p1[0], p1[1], p1[2]);
    GEO::vec3 q2(p2[0], p2[1], p2[2]);
    return sf_tree->segment_intersection(q1, q2);
  }

  // Finds the intersection between a segment and a surface that
  // is nearest to the first extremity (p1) of the segment.
  //
  // Note: if no intersection: store Vector3::Zero();
  inline bool sf_segment_nearest_intersection(const Vector3 &p1,
                                              const Vector3 &p2,
                                              Vector3 &intersection) const {
    intersection = Vector3::Zero();
    GEO::vec3 q1(p1[0], p1[1], p1[2]);
    GEO::vec3 q2(p2[0], p2[1], p2[2]);
    // [out]	t	if there was an intersection, it is t*q2 + (1-t)*q1
    // [out]	f	the intersected nearest facet or index_t(-1) if there
    // was no intersection.
    double t;
    GEO::index_t f;
    bool is_intersect = sf_tree->segment_nearest_intersection(q1, q2, t, f);
    if (is_intersect) {
      auto inter = t * q2 + (1 - t) * q1;
      intersection = Vector3(inter[0], inter[1], inter[2]);
    }
    return is_intersect;
  }

 public:
  // check if given polygon intersects sf_mesh
  std::set<int> sf_polygon_intersection(
      const GEO::Mesh &sf_mesh, const std::vector<Vector3> &polygon) const;
  // return all facets of sf_mesh whose bbox intersects the bbox of given
  // polygon (NOTE: not intersect with polygon itself, too slow to implement)
  std::vector<int> sf_polygon_box_intersections(
      const GEO::Mesh &sf_mesh, const std::vector<Vector3> &polygon) const;

 public:
  ///////////////////
  // for sharp edges
  inline double project_to_s(Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //?
    s_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];

    return sq_dist;
  }

  ///////////////////
  // for sf_mesh
  inline int project_to_sf_get_nearest_face(Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];
    return fidx;
  }

  inline Scalar project_to_sf(Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];

    return sq_dist;
  }
  inline int get_nearest_face_sf(const Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
  }

  inline int get_nearest_face_sf(const Vector3 &p, double &sq_dist) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    return sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
  }

  inline int get_nearest_point_on_sf(const Vector3 &p, Vector3 &q,
                                     double &sq_dist) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    int fidx = sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    q = Vector3(nearest_p[0], nearest_p[1], nearest_p[2]);
    return fidx;
  }

  inline Scalar get_sq_dist_to_sf(const Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    sf_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    return sq_dist;
  }

  inline Scalar get_sq_dist_to_s(const Vector3 &p) const {
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    s_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    return sq_dist;
  }

  /////////////////////
  // for concave lines
  inline int project_to_cc_get_nearest_face(Vector3 &p, double &sq_dist) const {
    if (!is_cc_mesh_exist) return -1;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    int eid = cc_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    p[0] = nearest_p[0];
    p[1] = nearest_p[1];
    p[2] = nearest_p[2];
    return eid;
  }

  inline int get_nearest_edge_cc(const Vector3 &p) const {
    if (!is_cc_mesh_exist) return -1;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    return cc_tree->nearest_facet(geo_p, nearest_p, sq_dist);
  }

  inline int get_nearest_edge_cc(const Vector3 &p, double &sq_dist) const {
    if (!is_cc_mesh_exist) return -1;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    sq_dist = std::numeric_limits<double>::max();  //??
    return cc_tree->nearest_facet(geo_p, nearest_p, sq_dist);
  }

  inline Scalar get_sq_dist_to_cc(const Vector3 &p) const {
    if (!is_cc_mesh_exist) return DBL_MAX;
    GEO::vec3 geo_p(p[0], p[1], p[2]);
    GEO::vec3 nearest_p;
    double sq_dist = std::numeric_limits<double>::max();  //??
    cc_tree->nearest_facet(geo_p, nearest_p, sq_dist);
    return sq_dist;
  }

  //// closeness check - point
  inline bool is_close_to_sf(const Vector3 &p, const Scalar eps_2) const {
    double sq_dist = get_sq_dist_to_sf(p);
    if (Scalar(sq_dist) < eps_2) return true;
    return false;
  }
  inline bool is_close_to_s(const Vector3 &p, const Scalar eps_2,
                            double &sq_dist) const {
    sq_dist = get_sq_dist_to_s(p);
    if (Scalar(sq_dist) < eps_2) return true;
    return false;
  }
};

}  // namespace matfp
