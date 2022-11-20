// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/MedialMeshGenerator.h"

#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>

#include "matfp/InscribedSpheres.h"
#include "matfp/MeshProcessor.h"
#include "matfp/RPDGenerator.h"
#include "matfp/Triangulation.h"
#include "matfp/WindingFilter.h"

#define USE_MVERTEX_FOR_MAT true

namespace matfp {

bool is_delete_mat_if_cover_sharp_edge(
    const std::vector<double>& seeds, const AABBWrapper& aabb_wrapper,
    const Vector3& center, const double& sq_radius,
    const std::set<std::array<int, 2>>& k_nearest_se, bool is_debug) {
  /*
   * http://paulbourke.net/geometry/circlesphere/index.html#linesphere
   * Calculate the intersection of a segment and a sphere
   * The line segment is defined from p1 to p2
   * The sphere is of radius r and centered at O
   * There are potentially point p of intersection given by
   * p = p1 + mu1 (p2 - p1)
   * also if p is the closest point on segment to sphere
   * (p-O).dot(p2-p1) = 0
   * Return TRUE if mu in [0, 1] and |Op| <= r.
   */
  auto is_one_se_covered_by_sphere =
      [&](const Vector3& p1, const Vector3& p2, const Vector3& O,
          const double& sq_radius, const bool& is_debug) {
        Vector3 p1p2 = p2 - p1;
        double mu = (O - p1).dot(p1p2) / p1p2.dot(p1p2);
        Vector3 p = p1 + mu * (p2 - p1);
        if (is_debug) logger().debug("mu: {}", mu);

        if (mu <= 1. && mu >= 0.) {
          double Op_len = (p - O).norm();
          double radius = std::sqrt(sq_radius);
          if (is_debug) logger().debug("Op_len: {}, radius {}", Op_len, radius);
          if (Op_len <= radius) {
            return true;
          }
        }
        return false;
      };

  // using aabb_wrapper is necessary
  // this would avoid many errors when adding new points to sharp edges
  double dist = std::sqrt(std::fabs(aabb_wrapper.get_sq_dist_to_s(center)));
  double delta = dist - std::sqrt(std::fabs(sq_radius));
  if (delta < SCALAR_ZERO) {
    if (is_debug)
      logger().debug("aabb_wrapper delta {} < SCALAR_ZERO, delete", delta);
    return true;
  } else {
    bool is_delete = false;
    for (const auto one_edge : k_nearest_se) {
      int se_idx = one_edge[0];
      int adj_se_i = one_edge[1];
      Vector3 v0(seeds[se_idx * 3], seeds[se_idx * 3 + 1],
                 seeds[se_idx * 3 + 2]);
      Vector3 v1(seeds[adj_se_i * 3], seeds[adj_se_i * 3 + 1],
                 seeds[adj_se_i * 3 + 2]);
      is_delete =
          is_one_se_covered_by_sphere(v0, v1, center, sq_radius, is_debug);
      if (is_debug)
        logger().debug("one_edge: {}, is_delete: {}", one_edge, is_delete);
    }
    return is_delete;
  }
};

void create_mat_vertices_from_RT(const std::vector<MVertex>& all_medial_spheres,
                                 const std::vector<int>& valid_medial_spheres,
                                 const RegularTriangulationNN& rt,
                                 NonManifoldMesh& mat, bool is_debug) {
  if (is_debug) logger().debug("creating mat vertices ...");
  mat.vertices.clear();
  mat.numVertices = 0;  // updated in mat.create_vertex()
  // using valid_medial_spheres here because Finite_vertices_iterator_rt may
  // iterate vertices randomly, not by valid_idx
  for (int valid_idx = 0; valid_idx < valid_medial_spheres.size();
       valid_idx++) {
    const int all_idx = valid_medial_spheres[valid_idx];
    const MVertex& mat_p = all_medial_spheres.at(all_idx);
    if (mat_p.is_outside || mat_p.is_on_bbox) {
      logger().error("mat_p {} should not be MAT vertex", mat_p.tag);
      mat_p.print_info();
      log_and_throw("ERROR");
    };
    // always store the original radius
    double radius = std::sqrt(mat_p.sq_radius);
    mat.create_vertex(mat_p.pos[0], mat_p.pos[1], mat_p.pos[2], radius,
                      valid_idx, false /*is_outside*/, mat_p.type, all_idx);
    // sanity check
    if (valid_idx != mat.numVertices - 1) {
      logger().error(
          "ERROR creating mat vertex, valid_idx: {} != mat.numVertices-1: {}",
          valid_idx, mat.numVertices - 1);
      log_and_throw("ERROR");
    }
  }
}

void create_mat_tets_from_RT(const RegularTriangulationNN& rt,
                             NonManifoldMesh& mat, bool is_debug) {
  if (is_debug) logger().debug("creating mat tets if any ...");
  // if tet's dual point is inside, then this tet is suppose to be in MAT
  // (even though we do not want this)
  for (const RTI& rti : rt.rt_ts_info) {
    if (rti.is_dual_point_outside) continue;
    std::vector<int> vids;
    std::copy(rti.valid_vs_tags.begin(), rti.valid_vs_tags.end(),
              std::back_inserter(vids));
    // // skip creating tet if contains a corner
    // if (rti.is_contain_corner) {
    //   // logger().debug("skip creating tet for all_tags: {}",
    //   // rti.all_vs_tags);
    //   // continue;
    // }
    int tid = mat.create_tet(vids, -1 /*tag*/, true /*is_auto_create_face*/);
    if (tid == -1) log_and_throw("tet must be created!!");
    mat.tets[tid].dual_vertex = rti.center;
    mat.tets[tid].is_dual_vertex_outside = rti.is_dual_point_outside;
    mat.tets[tid].dual_sq_radius = rti.sq_radius;
    mat.tets[tid].all_vs_tags = rti.all_vs_tags;
  }  // for rt.rt_ts_info
}

void save_mat_face_dual_info(const RFI& rfi, NonManifoldMesh_Face& face) {
  face.rfi_id = rfi.id;
  face.all_vs_tags = rfi.all_vs_tags;
  if (rfi.all_vs_tags[0] == -1) {
    rfi.print_info();
    log_and_throw("ERROR: save_mat_face_dual_info");
  }
  face.dual_segment = rfi.dual_segment;
  face.is_dual_seg_vs_outside = rfi.is_dual_seg_vs_outside;
  face.is_dual_seg_intersect = rfi.is_dual_seg_intersect;
  face.dual_intersections = rfi.dual_intersections;
  face.dist_dual_intersections = rfi.dist_dual_intersections;
}

void save_mat_edge_dual_info(const REI& rei, NonManifoldMesh_Edge& edge) {
  edge.rei_id = rei.id;
  edge.all_vs_tags = rei.all_vs_tags;
  edge.dual_polygon = rei.dual_polygon;
  edge.is_dual_poly_vs_outside = rei.is_dual_poly_vs_outside;
  edge.is_dual_poly_intersect = rei.is_dual_poly_intersect;
  edge.poly_bbox_isect_triangles = rei.poly_bbox_isect_triangles;
}

void create_mat_faces_from_RT(const std::vector<MVertex>& all_medial_spheres,
                              const RegularTriangulationNN& rt,
                              NonManifoldMesh& mat, bool is_debug) {
  if (is_debug) logger().debug("creating mat faces ...");
  for (const RFI& rfi : rt.rt_fs_info) {
    // could be -1 (8 bbox spheres)
    const std::array<int, 3>& vvid = rfi.valid_vs_tags;
    const std::array<int, 3>& all_tags = rfi.all_vs_tags;
    const std::array<Vector3, 2>& f_dual_segment = rfi.dual_segment;

    // NOTE: this need to be checked first!!!
    // face already created (by tets)
    // we only store dual info now
    int fid;
    std::set<int> vvid_set;
    vvid_set.insert(vvid.begin(), vvid.end());
    if (mat.Face(vvid_set, fid)) {
      save_mat_face_dual_info(rfi, *(mat.faces[fid].second));
      continue;
    }
    // sanity check
    for (const auto& adj_tid : rfi.adj_tets) {
      if (rt.rt_ts_info[adj_tid].is_contain_corner) continue;  // skip checking
      if (!rt.rt_ts_info[adj_tid].is_dual_point_outside) {
        // face must be created already
        logger().error(
            "ERROR: mat face must be created already for rt_ts_info {} that is "
            "inside",
            adj_tid);
        log_and_throw("ERROR");
      }
    }

    // if (!rfi.is_dual_seg_vs_outside[0] || !rfi.is_dual_seg_vs_outside[1]) {
    //   logger().error("ERROR: rfi face must already been created in tet!");
    //   rfi.print_info();
    //   log_and_throw("ERROR");
    // }

    bool is_create_mat_face = true;
    // both adjacent tets are outside
    // check face's dual segment
    // if dual segment has no intersection with input mesh, then skip
    if (!rfi.is_dual_seg_intersect) is_create_mat_face = false;
    // for other reasons, maybe
    // 1. any vertex is outside
    // 2. overlap by other spheres
    // 3. maybe be bbox spheres
    for (int i = 0; i < rfi.all_vs_tags.size(); i++) {
      int all_idx = rfi.all_vs_tags[i];
      int valid_idx = rfi.valid_vs_tags[i];
      if (valid_idx < 0 || all_idx < 0) {
        is_create_mat_face = false;
        break;
      }
      const auto& mat_p = all_medial_spheres[all_idx];
      if (all_idx > 0 && mat_p.is_outside) {
        is_create_mat_face = false;
        break;
      }
    }
    if (!is_create_mat_face) continue;

    // if a feature sphere is not connecting to a neighboring feature sphere
    int cnt_features = 0;
    for (const auto& all_idx1 : all_tags) {
      const MVertex& mat_p1 = all_medial_spheres.at(all_idx1);
      if (!mat_p1.is_a_feature_sphere()) continue;
      cnt_features++;
      // for (const auto& all_idx2 : all_tags) {
      //   if (all_idx1 == all_idx2) continue;
      //   const MVertex& mat_p2 = all_medial_spheres.at(all_idx2);
      //   if (!mat_p2.is_a_feature_sphere()) continue;
      //   if (mat_p1.se_adj_se.find(all_idx2) == mat_p1.se_adj_se.end()) {
      //     // not neighboring on the same sharp edge
      //     is_create_mat_face = false;
      //     break;
      //   }
      // }
      // if (!is_create_mat_face) continue;
    }
    if (!is_create_mat_face) continue;
    // 3 vertices cannot all on feature!!!
    // if (cnt_features == 3) continue;

    // logger().debug("creating face: {}", vvid);
    fid = mat.create_face(vvid_set);
    save_mat_face_dual_info(rfi, *(mat.faces[fid].second));

    // logger().error("mat face {} has tags {} is_dual_seg_intersect {},
    // is_create_mat_face {}, dual segment ({},{},{}), ({},{},{})",
    //     fid, all_tags, rfi.is_dual_seg_intersect, is_create_mat_face,
    //     f_dual_segment[0][0], f_dual_segment[0][1], f_dual_segment[0][2],
    //     f_dual_segment[1][0], f_dual_segment[1][1], f_dual_segment[1][2]
    // );
  }  // for rt.rt_fs_info

  if (is_debug) logger().debug("creating mat faces around corners ...");
  // create corner MAT faces
  for (const auto& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted) continue;
    if (!mat_p.is_added_for_corner) continue;
    if (mat_p.valid_idx < 0) continue;
    // added_for_two_spheres could contains {-1,-1} for concave edges
    for (const auto& se_tag_pair : mat_p.added_for_two_spheres) {
      std::set<int> vvid_set;
      vvid_set.insert(mat_p.valid_idx);
      for (const auto& one_se_tag : se_tag_pair) {
        if (one_se_tag < 0 || one_se_tag >= all_medial_spheres.size()) continue;
        const auto& mat_p_se = all_medial_spheres.at(one_se_tag);
        if (mat_p_se.valid_idx < 0) break;
        vvid_set.insert(mat_p_se.valid_idx);
      }
      if (vvid_set.size() != 3) continue;
      // create new MAT face around corner
      int fid;
      if (mat.Face(vvid_set, fid)) continue;
      fid = mat.create_face(vvid_set);
      if (is_debug)
        logger().debug("creating mat face around corner: {}", vvid_set);
    }  // for mat.added_for_two_spheres
  }    // for all_medial_spheres
}

void create_mat_edges_from_RT(const std::vector<MVertex>& all_medial_spheres,
                              const RegularTriangulationNN& rt,
                              NonManifoldMesh& mat, bool is_debug) {
  if (is_debug) logger().debug("creating mat edges ...");
  for (const REI& rei : rt.rt_es_info) {
    // vvid maybe -1 (8 bbox spheres)
    const std::array<int, 2>& vvid = rei.valid_vs_tags;
    const std::array<int, 2>& all_tags = rei.all_vs_tags;
    const std::vector<Vector3>& e_dual_polygon = rei.dual_polygon;

    // NOTE: check this first!!
    // already created
    int eid = -1;
    if (mat.Edge(vvid[0], vvid[1], eid)) {
      mat.eid2rei[eid] = rei.id;
      save_mat_edge_dual_info(rei, *(mat.edges[eid].second));
      continue;
    }

    bool is_create_mat_edge = true;
    // if dual polygon has no intersection with input, then skip
    if (!rei.is_dual_poly_intersect) is_create_mat_edge = false;
    // for other reasons, maybe
    // 1. any vertex is outside
    // 2. overlap by other spheres
    // 3. maybe be bbox spheres (valid_idx == -1)
    for (int i = 0; i < rei.all_vs_tags.size(); i++) {
      int all_idx = rei.all_vs_tags[i];
      int valid_idx = rei.valid_vs_tags[i];
      if (valid_idx < 0 || all_idx < 0) {
        is_create_mat_edge = false;
        break;
      }
      const auto& mat_p = all_medial_spheres[all_idx];
      if (all_idx > 0 && mat_p.is_outside) {
        is_create_mat_edge = false;
        break;
      }
    }
    if (!is_create_mat_edge) continue;

    // create a mat edge
    eid = mat.create_edge(vvid[0], vvid[1]);
    mat.eid2rei[eid] = rei.id;
    save_mat_edge_dual_info(rei, *(mat.edges[eid].second));
    // logger().debug("mat edge {} with all_idx {} not exist, created a new
    // one", eid, all_tags);
  }  // for rt.rt_es_info
}

void create_mat_from_RT(const std::vector<MVertex>& all_medial_spheres,
                        const std::vector<int>& valid_medial_spheres,
                        const RegularTriangulationNN& rt,
                        const std::string& mat_name, NonManifoldMesh& mat,
                        bool is_using_dilated_radius, bool is_debug) {
  logger().debug("creating mat from RT ....");
  mat.clear();
  mat.mat_name = mat_name;
  create_mat_vertices_from_RT(all_medial_spheres, valid_medial_spheres, rt, mat,
                              is_debug);
  create_mat_tets_from_RT(rt, mat, is_debug);
  create_mat_faces_from_RT(all_medial_spheres, rt, mat, is_debug);
  // create_mat_edges_from_RT(all_medial_spheres, rt, mat, is_debug);

  // check if vertices are active
  for (auto& vertex_tmp : mat.vertices) {
    auto& vertex = *(vertex_tmp.second);
    if (vertex.edges_.empty() && vertex.faces_.empty()) {
      vertex.is_deleted = true;
      mat.numVertices_active--;
    }
  }

  if (is_debug) logger().debug("creating normals and centroids ...");
  mat.ComputeFacesCentroid();
  mat.ComputeFacesNormal();
  mat.compute_edges_cone();
  mat.compute_faces_simple_triangles();
  mat.print_info();
}

}  // namespace matfp