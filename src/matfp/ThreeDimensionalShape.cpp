// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/ThreeDimensionalShape.h"

#include <CGAL/bounding_box.h>

#include "matfp/Common.h"
#include "matfp/Logger.h"
#include "matfp/WindingFilter.h"

namespace matfp {

SpecialEdge::SpecialEdge(const int _id, const EdgeType& _t,
                         const std::array<int, 2>& _vs_ids) {
  id = _id;
  type = _t;
  vs_seed_ids = _vs_ids;
}

void SpecialEdge::print_info() const {
  logger().debug(
      "-- SpecialEdge {}: type: {}, cc_line_id: {}, mat_se_tags: {}, "
      "vs_seed_ids: {}, "
      "adj_ref_fs_pair: {}",
      id, static_cast<int>(type), cc_line_id, mat_se_tags, vs_seed_ids,
      adj_ref_fs_pair);
}

SpecialCorner::SpecialCorner(const int _id, const int _s_id,
                             const Vector3& _pos) {
  id = _id;
  seed_id = _s_id;
  pos = _pos;
}

void SpecialCorner::print_info() const {
  logger().debug(
      "-- SpecialCorner {}: seed_id: {}, mat_tag: {}, adj_edges: {}, "
      "adj_edges_sorted: {}, "
      "pos: ({},{},{})",
      id, seed_id, mat_tag, adj_edges, adj_edges_sorted, pos[0], pos[1],
      pos[2]);
}

void SpecialCorner::push_to_adj_edges(const int sp_eid) {
  adj_edges.insert(sp_eid);
}

void MeshWrapper::get_fs_pair_given_edge(const std::array<int, 2>& e,
                                         std::vector<int>& fs_pair) const {
  fs_pair.clear();
  set_intersection(conn_tris.at(e.at(0)), conn_tris.at(e.at(1)), fs_pair);
  if (fs_pair.size() != 2) {
    logger().error("edge {} cannot find 2 adjacent faces: {}", e, fs_pair);
    log_and_throw("ERROR");
  }
}

Vector3 MeshWrapper::get_v_pos(const int vid) const {
  if (vid < 0 || vid >= input_vertices.size()) {
    logger().error("vid {} is out of range of sf_mesh [{},{})", vid, 0,
                   input_vertices.size());
    log_and_throw("Error");
  }
  return input_vertices[vid];
}

Vector3 MeshWrapper::get_f_normal(const int fid) const {
  if (fid < 0 || fid >= input_fnormals.size()) {
    logger().error("fid {} is out of range of sf_mesh [{},{})", fid, 0,
                   input_fnormals.size());
    log_and_throw("Error");
  }
  return input_fnormals[fid];
}

Vector3 MeshWrapper::get_f_centroid(const int fid) const {
  if (fid < 0 || fid >= input_faces.size()) {
    logger().error("fid {} is out of range of sf_mesh [{},{})", fid, 0,
                   input_faces.size());
    log_and_throw("Error");
  }
  return get_triangle_centroid(input_vertices[input_faces[fid][0]],
                               input_vertices[input_faces[fid][1]],
                               input_vertices[input_faces[fid][2]]);
}

void ThreeDimensionalShape::reload_s_edges_adjs() {
  s_edges_adjs.clear();
  // save s_edges in another way
  for (auto const& se_pair : s_edges) {
    s_edges_adjs[se_pair[0]].insert(se_pair[1]);
    s_edges_adjs[se_pair[1]].insert(se_pair[0]);
  }
}

void ThreeDimensionalShape::reload_ref_fs_pairs_not_cross() {
  // while calling iterate_sphere() to update tangent plane
  // fid cannot cross sharp edges or concave edges
  ref_fs_pairs_not_cross.clear();
  ref_fs_pairs_not_cross = se_ref_fs_pairs;
  for (const auto& one_cc_line : tan_cc_lines) {
    ref_fs_pairs_not_cross.insert(one_cc_line.adj_ref_fs_pair);
  }
}

}  // namespace matfp