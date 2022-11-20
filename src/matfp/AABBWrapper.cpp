// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/AABBWrapper.h"

#include <geogram/mesh/mesh_reorder.h>

#include "matfp/Logger.h"
#include "matfp/MeshProcessor.h"

namespace matfp {

bool AABBWrapper::init_mesh_from_edges(
    const std::vector<Vector3> &input_vertices,
    const std::vector<std::array<int, 2>> &edges, GEO::Mesh &mesh) {
  if (edges.empty()) {
    mesh.vertices.clear();
    mesh.vertices.create_vertices(1);
    mesh.vertices.point(0) = GEO::vec3(0, 0, 0);
    mesh.facets.clear();
    mesh.facets.create_triangles(1);
    mesh.facets.set_vertex(0, 0, 0);
    mesh.facets.set_vertex(0, 1, 0);
    mesh.facets.set_vertex(0, 2, 0);
    return false;
  } else {
    mesh.vertices.clear();
    mesh.vertices.create_vertices((int)edges.size() * 2);
    int cnt = 0;
    for (auto &e : edges) {
      for (int j = 0; j < 2; j++) {
        GEO::vec3 &p = mesh.vertices.point(cnt++);
        p[0] = input_vertices[e[j]][0];
        p[1] = input_vertices[e[j]][1];
        p[2] = input_vertices[e[j]][2];
      }
    }
    mesh.facets.clear();
    mesh.facets.create_triangles(
        (int)edges.size());  // degenerated facets -> edge
    for (int i = 0; i < edges.size(); i++) {
      mesh.facets.set_vertex(i, 0, i * 2);
      mesh.facets.set_vertex(i, 1, i * 2);
      mesh.facets.set_vertex(i, 2, i * 2 + 1);
    }
  }
  // mesh_reorder(mesh, GEO::MESH_ORDER_MORTON);
  return true;
}

void AABBWrapper::init_feature_meshes_and_trees(
    const std::vector<Vector3> &input_vertices,
    const std::vector<std::array<int, 2>> &s_edges,
    const std::vector<std::array<int, 2>> &cc_edges) {
  // for sharp edges
  s_mesh.clear(false, false);
  is_s_mesh_exist = init_mesh_from_edges(input_vertices, s_edges, s_mesh);
  s_tree = std::make_shared<GEO::MeshFacetsAABB>(s_mesh, true /*reorder*/);

  // for concave edges
  cc_mesh.clear(false, false);
  is_cc_mesh_exist = init_mesh_from_edges(input_vertices, cc_edges, cc_mesh);
  cc_tree = std::make_shared<GEO::MeshFacetsAABB>(cc_mesh, true /*reorder*/);
}

std::vector<int> AABBWrapper::sf_polygon_box_intersections(
    const GEO::Mesh &sf_mesh, const std::vector<Vector3> &polygon) const {
  // box of the polygon
  GEO::Box box;
  box.xyz_min[0] = box.xyz_max[0] = polygon[0][0];
  box.xyz_min[1] = box.xyz_max[1] = polygon[0][1];
  box.xyz_min[2] = box.xyz_max[2] = polygon[0][2];
  for (size_t j = 0; j < polygon.size(); j++) {
    for (int i = 0; i < 3; i++) {
      box.xyz_min[i] = std::min(box.xyz_min[i], polygon[j][i]);
      box.xyz_max[i] = std::max(box.xyz_max[i], polygon[j][i]);
    }
  }
  std::vector<int> poly_bbox_isect_triangles;
  auto action_tri_intersect = [&poly_bbox_isect_triangles](GEO::index_t f) {
    poly_bbox_isect_triangles.push_back(f);
  };

  sf_tree->compute_bbox_facet_bbox_intersections(box, action_tri_intersect);
  return poly_bbox_isect_triangles;
}

}  // namespace matfp