// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "sharp_feature_detection.h"

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/unique_rows.h>
#include <igl/writeOFF.h>
#include <matfp/Logger.h>

#include "pre_meshIO.h"

// to use logger
using namespace matfp;

namespace pre_matfp {

void mark_feature_attributes(const std::set<std::array<int, 2>>& s_edges,
                             const std::set<std::array<int, 2>>& cc_edges,
                             const std::set<int>& corners, GEO::Mesh& input) {
  logger().debug("Marking feature attributes ...");
  // Setup corners
  GEO::Attribute<int> attr_corners(input.vertices.attributes(), "corner");
  for (int i = 0; i < input.vertices.nb(); i++) {
    if (corners.find(i) != corners.end())
      attr_corners[i] = 1;
    else
      attr_corners[i] = 0;
  }

  // Setup sharp edges and concave edges
  GEO::Attribute<int> attr_se(input.edges.attributes(), "se");
  GEO::Attribute<int> attr_cce(input.edges.attributes(), "cce");
  for (int e = 0; e < input.edges.nb(); e++) {
    std::array<int, 2> edge = {
        {(int)input.edges.vertex(e, 0), (int)input.edges.vertex(e, 1)}};
    std::sort(edge.begin(), edge.end());
    if (s_edges.find(edge) != s_edges.end())
      attr_se[e] = 1;
    else if (cc_edges.find(edge) != cc_edges.end())
      attr_cce[e] = 1;
    else {
      attr_se[e] = 0;
      attr_cce[e] = 0;
    }
  }
}

// load GEO::Mesh from path
// reorder and mark feature attributes
bool load_mesh_and_preprocess(const Args& args, GEO::Mesh& input) {
  std::string path = args.input_surface_path;
  logger().debug("Loading mesh at {}...", path);
  input.clear(false, false);
  const bool ok = GEO::mesh_load(path, input);

  if (!ok) return false;
  if (!input.facets.are_simplices()) {
    mesh_repair(input, GEO::MeshRepairMode(GEO::MESH_REPAIR_TRIANGULATE |
                                           GEO::MESH_REPAIR_QUIET));
  }
  // This is important!
  // GEO::MeshFacetsAABB::initialize() requires mesh reorder
  // to make searching fast
  GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);

  std::vector<Vector3> points;
  std::vector<Vector3i> faces;
  normalize_mesh(input, points, faces);
  remove_duplicates(points, faces);
  if (args.num_subdivide != 0) {
    mesh_subdivision(points, faces, args.num_subdivide);
  }
  // reload and reorder both GEO::Mesh and points&facets
  reorder_mesh_from_vector(points, faces, input);

  // find feature edges and corners
  std::set<std::array<int, 2>> s_edges;
  std::set<std::array<int, 2>> cc_edges;
  std::set<int> corners;
  find_feature_edges(args, points, faces, s_edges, cc_edges, corners);
  mark_feature_attributes(s_edges, cc_edges, corners, input);

  logger().debug("#v: {}, #e: {}, #f: {}", input.vertices.nb(),
                 input.edges.nb(), input.facets.nb());

  return ok;
}

void find_feature_edges(const Args& args,
                        const std::vector<Vector3>& input_vertices,
                        const std::vector<Vector3i>& input_faces,
                        std::set<std::array<int, 2>>& s_edges,
                        std::set<std::array<int, 2>>& cc_edges,
                        std::set<int>& corners) {
  logger().debug(
      "Detecting sharp/concave edges and corners using threshold: {}",
      args.thres_concave);
  s_edges.clear();
  cc_edges.clear();
  corners.clear();

  std::vector<std::array<int, 2>> edges;
  std::map<int, std::unordered_set<int>> conn_tris;
  for (int i = 0; i < input_faces.size(); i++) {
    const auto& f = input_faces[i];
    for (int j = 0; j < 3; j++) {
      std::array<int, 2> e = {{f[j], f[(j + 1) % 3]}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
      conn_tris[input_faces[i][j]].insert(i);
    }
  }
  vector_unique(edges);

  // find sharp edges and concave edges
  for (const auto& e : edges) {
    std::vector<int> n12_f_ids;
    set_intersection(conn_tris[e[0]], conn_tris[e[1]], n12_f_ids);

    if (n12_f_ids.size() == 1) {  // open boundary
      logger().error("Detect open boundary!!! edge {} has only 1 face: {}", e,
                     n12_f_ids);
      log_and_throw("ERROR: we don't know how to handle open boundary!!");
    }
    int f_id = n12_f_ids[0];
    int j = 0;
    for (; j < 3; j++) {
      if ((input_faces[f_id][j] == e[0] &&
           input_faces[f_id][mod3(j + 1)] == e[1]) ||
          (input_faces[f_id][j] == e[1] &&
           input_faces[f_id][mod3(j + 1)] == e[0]))
        break;
    }
    Vector3 n = get_normal(input_vertices[input_faces[f_id][0]],
                           input_vertices[input_faces[f_id][1]],
                           input_vertices[input_faces[f_id][2]]);
    Vector3 c_n = get_triangle_centroid(input_vertices[input_faces[f_id][0]],
                                        input_vertices[input_faces[f_id][1]],
                                        input_vertices[input_faces[f_id][2]]);

    for (int k = 0; k < n12_f_ids.size(); k++) {
      if (n12_f_ids[k] == f_id) continue;
      Vector3 n1 = get_normal(input_vertices[input_faces[n12_f_ids[k]][0]],
                              input_vertices[input_faces[n12_f_ids[k]][1]],
                              input_vertices[input_faces[n12_f_ids[k]][2]]);
      Vector3 c_n1 =
          get_triangle_centroid(input_vertices[input_faces[n12_f_ids[k]][0]],
                                input_vertices[input_faces[n12_f_ids[k]][1]],
                                input_vertices[input_faces[n12_f_ids[k]][2]]);

      std::array<int, 2> ref_fs_pair = {{f_id, n12_f_ids[k]}};
      std::sort(ref_fs_pair.begin(), ref_fs_pair.end());
      std::array<Vector3, 2> ref_fs_normals = {{n, n1}};
      bool is_debug = false;

      //////////
      // Since cosine can only measure dihedral angle from (0, 180)
      // but concave has angle larger than 180
      // therefore we use different measurement for concave, and sharp edges
      //////////
      // Concave edges
      // c_n is a random vertex on plane A, c_n1 is a random vertex on plane B
      // n is normal of A
      //
      // 2021-09-04 ninwang:
      // If na and nb are the normals of the both adjacent faces,
      // and pa and pb vertices of the both faces that are not connected to
      // the edge, wherein na and pa belongs to the face A, and nb and pb to
      // the face B, then ( pb - pa ) . na > 0 => concave edge
      double tmp_concave = (c_n1 - c_n).normalized().dot(n);     // A, B
      double tmp_concave_2 = (c_n - c_n1).normalized().dot(n1);  // B, A
      if (tmp_concave > args.thres_concave ||
          tmp_concave_2 > args.thres_concave) {  // SCALAR_ZERO is too small
        if (is_debug)
          logger().debug("edge e {} is a concave edge, tmp_concave: {}", e,
                         tmp_concave);
        cc_edges.insert(e);  // once concave, never sharp
      } else {
        // Sharp edges (when it's convex)
        // angle between two normals of convex faces: theta
        // => cos(theta) = n1.dot(n)
        // acosine() range in [0, pi]
        // sharp edges => theta in (angle_sharp, pi)
        // here angle_sharp = 30
        // Note that, using theta CANNOT differentiate concave or convex,
        // so the concave detection must run first
        double tmp_convex = std::acos(n1.dot(n));
        double angle_sharp = PI * (args.thres_convex / 180.);
        if (angle_sharp < tmp_convex && tmp_convex < PI) {
          // logger().debug("sharp edge: theta is {}", tmp_convex);
          s_edges.insert(e);
        }
      }
    }  // for n12_f_ids
  }    // for edges

  // vector_unique(s_edges);

  // find corners
  // connect to at least 3 sharp edges
  // (not including concave edges)
  std::map<int, std::set<int>> neighbor_v;
  for (const auto& e : s_edges) {
    neighbor_v[e[0]].insert(e[1]);
    neighbor_v[e[1]].insert(e[0]);
  }
  for (const auto& pair : neighbor_v) {
    // Found a corner
    if (pair.second.size() > 2) {
      corners.insert(pair.first);
    }
  }

  logger().debug("#concave_edges = {}", cc_edges.size());
  logger().debug("#sharp_edges = {}", s_edges.size());
  logger().debug("#corners  = {}", corners.size());
}
}  // namespace pre_matfp