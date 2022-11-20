// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "Thinning.h"

#include "matfp/Common.h"

namespace matfp {

void add_face_edge_pairs(
    const NonManifoldMesh& mat, const std::set<int>& edges,
    const unsigned tet_id,  // which tet these edges belongs to
    std::queue<simple_pair>& queue, bool is_debug = false) {
  // find edges on only 1 boundary face
  for (const auto& eid : edges) {
    auto edge = *(mat.edges[eid].second);
    if (edge.is_deleted) continue;
    // if this edge is a mat feature, then skip
    int v1 = edge.vertices_.first;
    int v2 = edge.vertices_.second;
    if (mat.vertices[v1].second->type == 1 &&
        mat.vertices[v2].second->type == 1) {
      // logger().debug("[Prune] tid {} skip checking feature edge {} with
      // valid_idx ({},{})", tet_id, edge.tag, v1, v2);
      continue;
    }
    if (edge.faces_.size() == 1) {
      // push <face, edge> to queue
      queue.push(simple_pair(1, *(edge.faces_.begin()), eid, tet_id));
      if (is_debug) {
        logger().debug("[Prune] add tid {}, type {}: face {} edge {} to queue",
                       tet_id, queue.back().type, queue.back().idx0,
                       queue.back().idx1);
      }
    }
  }  // for mat edge
}

// just store importance for each mat face
void load_all_mat_face_importance_globally(NonManifoldMesh& mat,
                                           bool is_debug) {
  // for (const auto& one_tet : mat.tets) {
  //   bool is_debug = true;
  //   int tet_id = one_tet.tag;
  //   for (const auto& fid : one_tet.faces_) {
  //     auto& face = *(mat.faces[fid].second);
  for (const auto& face_tmp : mat.faces) {
    auto& face = *(face_tmp.second);
    // if (face.tets_.empty()) {
    //   // face.importance = DBL_MAX;
    //   log_and_throw("face of tet cannot have empty tets_!!!");
    //   continue;
    // }
    // do not count feature spheres
    // they are making avg_radius smaller
    double avg_radius = 0;
    int nb_extf = 0;
    for (const auto& vid : face.vertices_) {
      auto& vertex = *(mat.vertices[vid].second);
      if (vertex.type == SphereType::T_1_N) nb_extf++;
      avg_radius += vertex.radius;
    }
    int cnt = face.vertices_.size() - nb_extf;
    avg_radius /= std::max(cnt, 1);
    if (nb_extf > 2 || avg_radius <= SCALAR_ZERO) {
      // three vertices could all be feature vertices
      // this happends mostly for non-feature models
      face.importance = 0.01;
    } else {
      face.importance = face.dist_dual_intersections / (avg_radius * 2);
    }
  }  // for mat face
  // }    // for mat tets
  logger().debug("save mat faces importances");
}

// priority queue is ordered by first element of the pair
void sort_mat_face_importance_globally(NonManifoldMesh& mat,
                                       std::set<face_importance>& imp_queue,
                                       bool is_sort_randomly, bool is_debug) {
  imp_queue.clear();
  for (const auto& one_tet : mat.tets) {
    bool is_debug = true;
    int tet_id = one_tet.tag;
    for (const auto& fid : one_tet.faces_) {
      auto& face = *(mat.faces[fid].second);
      if (face.tets_.empty()) {
        // face.importance = DBL_MAX;
        log_and_throw("face of tet cannot have empty tets_!!!");
        continue;
      }
      if (is_sort_randomly) {
        face.importance = (double)std::rand() / (RAND_MAX);
        imp_queue.insert(std::make_pair(face.importance, fid));
      } else {
        // do not count feature spheres
        // they are making avg_radius smaller
        double avg_radius = 0;
        int nb_extf = 0;
        for (const auto& vid : face.vertices_) {
          auto& vertex = *(mat.vertices[vid].second);
          if (vertex.type == SphereType::T_1_N) nb_extf++;
          avg_radius += vertex.radius;
        }
        int cnt = face.vertices_.size() - nb_extf;
        avg_radius /= std::max(cnt, 1);
        if (nb_extf > 2 || avg_radius <= SCALAR_ZERO) {
          // three vertices could all be feature vertices
          // this happends mostly for non-feature models
          face.importance = 0.01;
        } else {
          face.importance = face.dist_dual_intersections / (avg_radius * 2);
        }
        imp_queue.insert(std::make_pair(face.importance, fid));
        // if (is_debug) {
        //   logger().debug(
        //       "[Prune] fid {} has importance: {}, pushed to priority queue",
        //       fid, face.importance);
        // }
      }
    }  // for mat face
  }    // for mat tets
  if (is_debug)
    logger().debug("[Prune] total {} mat faces in priority queue",
                   imp_queue.size());
}

void prune_one_simple_pair(const simple_pair& sp_del, NonManifoldMesh& mat,
                           bool is_debug) {
  if (is_debug) {
    logger().debug("[Prune] prune simple pair:");
    sp_del.print_info();
  }
  if (sp_del.type == 2) {  // tet-face pair
    int tet_id = sp_del.idx0;
    int fid = sp_del.idx1;
    mat.delete_tet(tet_id);
    mat.delete_face(fid);
    // if (is_debug) {
    //   auto& face = *(mat.faces[fid].second);
    //   face.delete_for_tet = tet_id;
    //   logger().debug("[Prune] deleted tet {} and face {}", tet_id, fid);
    //   logger().debug("[Prune] face {} have tets: {}, tets_old: {}", fid,
    //                  face.tets_, face.tets_old_);
    // }
  } else if (sp_del.type == 1) {  // face-edge pair
    int fid = sp_del.idx0;
    int eid = sp_del.idx1;
    mat.delete_face(fid);
    mat.delete_edge(eid);
    if (is_debug)
      logger().debug("[Prune] deleted face {} and edge {}", fid, eid);
  } else {
    // error
    logger().error("[Prune] unkown type {} to prune", sp_del.type);
    log_and_throw("ERROR");
  }
}

void prune_tets_while_iteration(NonManifoldMesh& mat,
                                std::set<face_importance>& imp_queue,
                                bool is_debug) {
  if (mat.tets.empty()) return;
  if (imp_queue.empty()) {
    logger().error("[Prune] imp_queue is not empty?? imp_queue {}",
                   imp_queue.size());
    return;
  }

  while (true) {
    if (is_debug)
      logger().debug("still {} mat tets to prune", mat.numTets_active);
    for (const auto& imp_pair : imp_queue) {
      double f_imp = imp_pair.first;
      int fid = imp_pair.second;
      auto& face = *(mat.faces[fid].second);
      if (face.is_deleted) continue;
      if (face.tets_.empty() || face.tets_.size() > 1) continue;
      if (face.tets_.size() != 1)
        log_and_throw("ERROR face has more than 1 tets");
      // everytime we delete a face, we break for loop of imp_queue
      // and checking again from face with smallest importance
      int tid = *(face.tets_.begin());
      // push <tet, face> to delete
      auto sp_del = simple_pair(2, tid, fid, tid);
      if (is_debug) {
        logger().debug(
            "[Prune] add tet_face pair to prune: tid {}, fid: {}, importance: "
            "{}",
            tid, fid, f_imp);
      }
      prune_one_simple_pair(sp_del, mat, is_debug);
      break;
    }  // for imp_queue

    // break while loop once we have all tets cleaned
    if (mat.numTets_active == 0) break;
  }  // while true
}

void prune_faces_while_iteration(const std::vector<MVertex>& all_medial_spheres,
                                 NonManifoldMesh& mat,
                                 std::set<face_importance>& imp_queue,
                                 double imp_thres, bool is_debug) {
  if (imp_queue.empty()) {
    logger().error("[Prune] imp_queue is not empty?? imp_queue {}",
                   imp_queue.size());
    return;
  }
  auto is_edge_on_ext_feature = [&](const std::pair<int, int>& vs_pair) {
    int all_idx1 = mat.vertices[vs_pair.first].second->all_idx;
    int all_idx2 = mat.vertices[vs_pair.second].second->all_idx;
    const auto mat_p1 = all_medial_spheres.at(all_idx1);
    const auto mat_p2 = all_medial_spheres.at(all_idx2);
    if (mat_p1.is_a_feature_sphere() && mat_p2.is_a_feature_sphere() &&
        mat_p1.se_adj_se.find(all_idx2) != mat_p1.se_adj_se.end() &&
        mat_p2.se_adj_se.find(all_idx1) != mat_p2.se_adj_se.end()) {
      return true;
    }
    // if ((mat.vertices[vs_pair.first].second)->type == SphereType::T_1_N &&
    //     (mat.vertices[vs_pair.second].second)->type == SphereType::T_1_N) {
    //   return true;
    // }
    return false;
  };

  int num_tet_faces_are_good = 0;
  while (true) {
    if (is_debug)
      logger().debug("num_tet_faces_are_good: {}, imp_queue.size {}",
                     num_tet_faces_are_good, imp_queue.size());
    num_tet_faces_are_good = 0;
    for (const auto& imp_pair : imp_queue) {
      bool is_break_queue_loop = false;
      double f_imp = imp_pair.first;
      int fid = imp_pair.second;
      auto& face = *(mat.faces[fid].second);
      if (face.importance != f_imp) {
        logger().error(
            "ERROR: face {} cannot have inconsistant importace: {} and {}",
            face.importance, f_imp);
        log_and_throw("ERROR");
      }
      // if (face.is_deleted) {
      if (face.is_deleted || face.importance >= imp_thres) {
        num_tet_faces_are_good++;
        continue;
      }
      if (!face.tets_.empty())
        log_and_throw("Tets of mat face cannot be empty");
      // everytime we delete a face, we break for loop of imp_queue
      // and checking again from face with smallest importance
      for (const auto& eid : face.edges_) {
        const auto edge = *(mat.edges[eid].second);
        if (edge.faces_.size() > 1) continue;  // check next edge
        // if edge is on sharp edge, then skip
        if (is_edge_on_ext_feature(edge.vertices_)) continue;
        // push <face, edge> to delete
        auto sp_del = simple_pair(1, fid, eid, -1);
        if (is_debug) {
          logger().debug(
              "[Prune] add face_edge pair to prune: fid {}, eid: {}, "
              "importance: "
              "{}",
              fid, eid, f_imp);
        }
        prune_one_simple_pair(sp_del, mat, is_debug);
        is_break_queue_loop = true;
        break;
      }  // for face.edges_
      if (is_break_queue_loop) break;
      // the face is good, every edge has 2 neighboring faces
      num_tet_faces_are_good++;
    }  // for imp_queue

    // break while loop once we have all tets faces cleaned
    if (num_tet_faces_are_good == imp_queue.size()) break;
  }  // while true

  return;
}

// imp_thres: more faces will be prunes if bigger
void prune(const std::vector<int>& valid_medial_spheres,
           const std::vector<MVertex>& all_medial_spheres, NonManifoldMesh& mat,
           double imp_thres, bool is_sort_randomly, bool is_debug) {
  logger().debug("start thinning with importance threshold: {} ...", imp_thres);
  // using set as priority queue
  // (importance, fid) starts with the smallest importance
  // tie-broken with smallest fid.
  std::set<face_importance> imp_queue;  // in ascending order
  sort_mat_face_importance_globally(mat, imp_queue, is_sort_randomly, is_debug);
  prune_tets_while_iteration(mat, imp_queue, is_debug);
  prune_faces_while_iteration(all_medial_spheres, mat, imp_queue, imp_thres,
                              is_debug);
}

}  // namespace matfp