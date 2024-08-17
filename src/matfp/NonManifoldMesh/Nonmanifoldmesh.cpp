// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "Nonmanifoldmesh.h"

#include <igl/writeOFF.h>
#include <igl/writePLY.h>

#include <boost/lexical_cast.hpp>
#include <cstdio>
#include <ctime>
#include <string>

#include "matfp/Common.h"

// ninwang: 2019-08-24
// this is for eliminate error BOOST_PARAMETER_MAX_ARITY must be at least 12 for
// CGAL::Mesh_3" while compiling
#define BOOST_PARAMETER_MAX_ARITY 20

namespace matfp {

void NonManifoldMesh_Vertex::print_mat_vertex() {
  logger().debug("------ Medial Vertex {} Info ------", tag);
  logger().debug("type: {}, is_deleted: {}", (int)type, is_deleted);
  logger().debug("pos: ({},{},{})", pos[0], pos[1], pos[2]);
  logger().debug("radius: {}", radius);
}

void NonManifoldMesh_Face::print_mat_face() const {
  logger().debug("------ MAT face {}, rfi_id: {}", tag, rfi_id);
  logger().debug("is_deleted {}, all_vs_tags: {}, is_dual_seg_intersect: {}",
                 is_deleted, all_vs_tags, is_dual_seg_intersect);
  logger().debug(
      "segment: ({},{},{}) -> ({},{},{}), is_dual_seg_vs_outside: {} ",
      dual_segment[0][0], dual_segment[0][1], dual_segment[0][2],
      dual_segment[1][0], dual_segment[1][1], dual_segment[1][2],
      is_dual_seg_vs_outside);
  logger().debug(
      "dual_intersection: ({},{},{}) -> ({},{},{}), dist_dual_intersections: "
      "{}, importance: {}",
      dual_intersections[0][0], dual_intersections[0][1],
      dual_intersections[0][2], dual_intersections[1][0],
      dual_intersections[1][1], dual_intersections[1][2],
      dist_dual_intersections, importance);
  logger().debug(
      "tets: {}, tets_old_: {}, delete_for_tet: {}, vertices {}, edges {}",
      tets_, tets_old_, delete_for_tet, vertices_, edges_);
}

void NonManifoldMesh_Edge::print_mat_edge() const {
  logger().debug("------ MAT edge {}, rei_id: {}", tag, rei_id);
  logger().debug("is_deleted {}, all_vs_tags: {}", is_deleted, all_vs_tags);
  logger().debug("is_dual_poly_vs_outside: {}, is_dual_poly_intersect: {}",
                 is_dual_poly_vs_outside, is_dual_poly_intersect);
  logger().debug("poly_bbox_isect_triangles: {}", poly_bbox_isect_triangles);
}

void NonManifoldMesh_Tet::print_mat_tet() const {
  logger().debug("------ MAT tet {}:", tag);
  logger().debug("all_vs_tags: {}, neigh_tets_: {}", all_vs_tags, neigh_tets_);
  logger().debug(
      "dual vertex: ({},{},{}), is_dual_vertex_outside: {}, dual_sq_radius: {}",
      dual_vertex[0], dual_vertex[1], dual_vertex[2], is_dual_vertex_outside,
      dual_sq_radius);
}

void NonManifoldMesh::print_info() {
  logger().info("------ Medial Mesh Info ------");
  logger().info("mat numVertices: {}, numVertices_active {}", numVertices,
                numVertices_active);
  logger().info("mat numEdges: {}, numEdges_active {}", numEdges,
                numEdges_active);
  logger().info("mat numFaces: {}, numFaces_active {}", numFaces,
                numFaces_active);
  logger().info("mat numTets: {}, numTets_active {}", numTets, numTets_active);
}

bool NonManifoldMesh::is_empty() {
  return vertices.empty() && edges.empty() && faces.empty();
}

int NonManifoldMesh::get_sphere_all_idx(const int& valid_idx) {
  if (valid_idx < 0 || valid_idx > vertices.size() - 1) {
    logger().debug("ERROR: valid_idx {} not in range [{},{}]", valid_idx, 0,
                   vertices.size() - 1);
    log_and_throw("ERROR");
  }
  return vertices[valid_idx].second->all_idx;
}

int NonManifoldMesh::create_vertex(const double x, const double y,
                                   const double z, const double radius,
                                   int tag, /*if not given then auto-increase*/
                                   const bool is_outside, const SphereType type,
                                   const int all_idx) {
  if (tag == -1) tag = numVertices;
  Bool_VertexPointer bvp;
  bvp.first = true;
  bvp.second = new NonManifoldMesh_Vertex;
  (*bvp.second).tag = tag;  // give each mat point an index = valid_idx
  (*bvp.second).is_outside = is_outside;
  (*bvp.second).pos[0] = x;
  (*bvp.second).pos[1] = y;
  (*bvp.second).pos[2] = z;
  (*bvp.second).radius = radius;  // feature mat point
  (*bvp.second).type = type;
  (*bvp.second).all_idx = all_idx;
  vertices.push_back(bvp);
  numVertices++;
  numVertices_active++;
  max_radius = std::max(max_radius, radius);
  return tag;
}

int NonManifoldMesh::create_edge(const int vid0, const int vid1,
                                 int tag  // if not given then auto-increase
) {
  // edge already exists
  if (Edge(vid0, vid1, tag)) return tag;
  if (tag == -1) tag = numEdges;
  Bool_EdgePointer bep;
  bep.first = true;
  bep.second = new NonManifoldMesh_Edge;
  (*bep.second).tag = tag;  // give each mat point an index
  (*bep.second).vertices_.first = vid0;
  (*bep.second).vertices_.second = vid1;
  (*vertices[vid0].second).edges_.insert(tag);
  (*vertices[vid1].second).edges_.insert(tag);
  edges.push_back(bep);
  numEdges++;
  numEdges_active++;
  return tag;
}

int NonManifoldMesh::create_face(const std::set<int>& vvid,
                                 int tag  // if not given then auto-increase
) {
  if (Face(vvid, tag)) return tag;
  if (tag == -1) tag = numFaces;
  Bool_FacePointer bfp;
  bfp.first = true;
  bfp.second = new NonManifoldMesh_Face;
  (*bfp.second).tag = tag;
  (*bfp.second).vertices_ = vvid;
  std::set<int> eid;
  for (const auto& v1 : vvid) {
    for (const auto& v2 : vvid) {
      if (v1 == v2) continue;
      int one_eid = -1;
      if (Edge(v1, v2, one_eid))
        (*bfp.second).edges_.insert(one_eid);
      else {
        one_eid = create_edge(v1, v2);
        (*bfp.second).edges_.insert(one_eid);
      }
      if (one_eid == -1) {
        logger().error("[MAT] face not contain an edge with -1 index, vvid: {}",
                       vvid);
        log_and_throw("ERROR");
      }
      eid.insert(one_eid);
    }
  }
  if (eid.size() != 3) {
    logger().error("[MAT] face not contain 3 edges, vvid: {}, edges {}", vvid,
                   eid);
    log_and_throw("ERROR");
  }
  for (const auto& v : vvid) {
    (*vertices[v].second).faces_.insert(tag);
  }
  for (const auto& e : eid) {
    (*edges[e].second).faces_.insert(tag);
  }
  faces.push_back(bfp);
  numFaces++;
  numFaces_active++;
  return tag;
}

int NonManifoldMesh::create_tet(const std::vector<int> vids, int tag,
                                bool is_auto_create_face) {
  if (vids.size() != 4) {
    logger().error("[MAT TET] fail to create MAT tet with vids: {}", vids);
    return -1;
  }
  if (tag == -1) tag = numTets;
  NonManifoldMesh_Tet tet;
  bool is_debug = false;
  tet.tag = tag;
  std::vector<int> fids;
  // save faces
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      for (int k = j + 1; k < 4; k++) {
        std::set<int> vvid;
        vvid.insert(vids[i]);
        vvid.insert(vids[j]);
        vvid.insert(vids[k]);
        int fid = -1;
        if (Face(vvid, fid)) {  // check if mat face exist
          // nothing
        } else {
          if (!is_auto_create_face) return -1;
          // creating new face if not exist
          fid = create_face(vvid);
        }
        const NonManifoldMesh_Face& face = *(faces[fid].second);
        tet.vertices_.insert(face.vertices_.begin(), face.vertices_.end());
        tet.edges_.insert(face.edges_.begin(), face.edges_.end());
        tet.faces_.insert(fid);
        tet.neigh_tets_.insert(face.tets_.begin(), face.tets_.end());
        fids.push_back(fid);
      }  // k
    }  // j
  }  // i

  if (fids.size() != 4) {
    logger().error("[MAT TET] error creating tet {}, fids size not 4: {}", tag,
                   fids);
    log_and_throw("ERROR");
  }
  if (is_debug) logger().debug("tet {} has faces {}", tag, fids);
  // only insert tid to face when tet is truly created.
  for (const auto& fid : fids) {
    NonManifoldMesh_Face& face = *(faces[fid].second);
    face.tets_.insert(tag);
    face.tets_old_.insert(tag);
  }
  tets.push_back(tet);
  numTets++;
  numTets_active++;
  return tag;
}

void NonManifoldMesh::clear() {
  for (auto& v_pointer : vertices) {
    delete v_pointer.second;
    v_pointer.second = NULL;
  }
  vertices.clear();

  for (auto& e_pointer : edges) {
    delete e_pointer.second;
    e_pointer.second = NULL;
  }
  edges.clear();

  for (auto& f_pointer : faces) {
    delete f_pointer.second;
    f_pointer.second = NULL;
  }
  faces.clear();

  // this will be cleaned later
  tets.clear();

  numVertices = numEdges = numFaces = numTets = 0;
  numVertices_active = numEdges_active = numFaces_active = numTets_active = 0;

  vs_adj_pair_to_skip.clear();

  // clear feature edges
  mat_intf_edges.clear();
  mat_intf_adjs.clear();
  mat_extf_edges.clear();

  // clear invalid edge
  invalid_mat_edges.clear();
}

int NonManifoldMesh::clear_all_tets() {
  tets_in_mat.clear();
  tets.clear();
  numTets = 0;
  numTets_active = 0;

  // all adjacent info will also be cleaned
  for (int fid = 0; fid < numFaces; fid++) {
    faces[fid].second->tets_.clear();
  }
}

bool NonManifoldMesh::delete_edge(const int eid) {
  if (eid < 0 || eid > numEdges) return false;
  auto& edge = *(edges[eid].second);
  if (edge.is_deleted) return true;
  edge.is_deleted = true;
  numEdges_active--;
  // update vertices
  vertices[edge.vertices_.first].second->edges_.erase(eid);
  if (vertices[edge.vertices_.first].second->edges_.empty())
    vertices[edge.vertices_.first].second->is_deleted = true;
  vertices[edge.vertices_.second].second->edges_.erase(eid);
  if (vertices[edge.vertices_.second].second->edges_.empty())
    vertices[edge.vertices_.second].second->is_deleted = true;
  // update faces
  for (const auto& fid : edge.faces_) {
    faces[fid].second->edges_.erase(eid);
  }
  return true;
}

bool NonManifoldMesh::delete_face(const int fid) {
  if (fid < 0 || fid > numFaces) return false;

  auto& face = *(faces[fid].second);
  if (face.is_deleted) return true;
  face.is_deleted = true;
  numFaces_active--;
  // update vertices
  for (const auto& vid : face.vertices_)
    vertices[vid].second->faces_.erase(fid);
  // update edges
  for (const auto& eid : face.edges_) {
    edges[eid].second->faces_.erase(fid);
    if (edges[eid].second->faces_.empty()) delete_edge(eid);
  }
  // update tets
  for (const auto& tid : face.tets_) tets[tid].faces_.erase(fid);
  return true;
}

bool NonManifoldMesh::delete_tet(const int tet_id) {
  if (tet_id < 0 || tet_id > numTets) return false;

  auto& tet = tets[tet_id];
  if (tet.is_deleted) return true;
  tet.is_deleted = true;
  numTets_active--;
  // update faces
  for (const auto& fid : tet.faces_) faces[fid].second->tets_.erase(tet_id);
  return true;
}

void NonManifoldMesh::get_adj_vertices(int vid,
                                       std::unordered_set<int>& adj_vs) const {
  adj_vs.clear();
  for (const auto& eid : vertices[vid].second->edges_) {
    const auto& mat_e = *(edges[eid].second);
    if (mat_e.is_deleted) continue;
    if (mat_e.vertices_.first != vid) adj_vs.insert(mat_e.vertices_.first);
    if (mat_e.vertices_.second != vid) adj_vs.insert(mat_e.vertices_.second);
  }
}

void NonManifoldMesh::add_vs_pair_to_skip(int v1, int v2) {
  // v1, v2 are mappint to valid_medial_spheres
  std::array<int, 2> vs_pair = {{v1, v2}};
  std::sort(vs_pair.begin(), vs_pair.end());
  vs_adj_pair_to_skip.insert(vs_pair);
}

bool NonManifoldMesh::is_skip_mat_edge_creation(int v1, int v2) {
  if (vs_adj_pair_to_skip.empty()) return false;

  std::array<int, 2> vs_pair = {{v1, v2}};
  std::sort(vs_pair.begin(), vs_pair.end());
  if (vs_adj_pair_to_skip.find(vs_pair) != vs_adj_pair_to_skip.end())
    return true;

  return false;
}

bool NonManifoldMesh::is_skip_mat_face_creation(const std::set<int>& vvid) {
  if (vs_adj_pair_to_skip.empty()) return false;

  for (const auto& v1 : vvid) {
    for (const auto& v2 : vvid) {
      if (v1 != v2 && is_skip_mat_edge_creation(v1, v2)) return true;
    }
  }

  // for (int i = 0; i < 3; i++) {
  // 	for (int j = i+1; j < 3; j++) {
  // 		if (is_skip_mat_edge_creation(vid[i], vid[j]))
  // 			return true;
  // 	}
  // }
  return false;
}

void NonManifoldMesh::computebb() {
  m_min[0] = 1e20;
  m_min[1] = 1e20;
  m_min[2] = 1e20;
  m_max[0] = -1e20;
  m_max[1] = -1e20;
  m_max[2] = -1e20;

  for (unsigned i = 0; i < vertices.size(); i++) {
    if (!vertices[i].first) continue;

    Vector3 ver = vertices[i].second->pos;
    if (ver[0] < m_min[0]) m_min[0] = ver[0];
    if (ver[1] < m_min[1]) m_min[1] = ver[1];
    if (ver[2] < m_min[2]) m_min[2] = ver[2];

    if (ver[0] > m_max[0]) m_max[0] = ver[0];
    if (ver[1] > m_max[1]) m_max[1] = ver[1];
    if (ver[2] > m_max[2]) m_max[2] = ver[2];
  }

  bbox_diag_l = std::sqrt(std::pow(m_max[0] - m_min[0], 2) +
                          std::pow(m_max[1] - m_min[1], 2) +
                          std::pow(m_max[2] - m_min[2], 2));
}

double NonManifoldMesh::get_nearest_distance(Vector3& input_p) {
  double min_dist = DBL_MAX;
  int min_index = -1;

  // loop MAT vertices
  for (int j = 0; j < numVertices; j++) {
    Vector3& mat_p = vertices[j].second->pos;
    double mat_r = vertices[j].second->radius;
    double tmp_length = std::abs((input_p - mat_p).norm() - mat_r);
    // update min_dist
    if (tmp_length < min_dist) {
      min_dist = tmp_length;
      min_index = j;  // store the index of mat vertex
    }
  }
  if (min_index == -1) {
    logger().debug(
        "ERROR: not found mat vertex for computing nearest distance");
    return min_dist;
  }
  // logger().debug("min_dist after vertices loop: {}", min_dist);

  // loop MAT edges
  set<int> near_edges = vertices[min_index].second->edges_;
  for (set<int>::iterator si = near_edges.begin(); si != near_edges.end();
       si++) {
    if (!edges[*si].first) continue;

    NonManifoldMesh_Edge mat_e = *edges[*si].second;

    NonManifoldMesh_Vertex v[2];
    Vector3 tfp;  // projection of input_p on edge {v0,v1}
    double td, tr;
    v[0] = *(vertices[mat_e.vertices_.first].second);
    v[1] = *(vertices[mat_e.vertices_.second].second);
    Vector3 v0 = v[0].pos;
    Vector3 v1 = v[1].pos;

    // check if input_p is within edge {v0,v1}, ignore if not
    double t((input_p - v0).dot(v1 - v0) / (v1 - v0).squaredNorm());
    if ((t >= 0.0) && (t <= 1.0)) {
      // projection of input_p on edge {v0,v1}
      tfp = (1.0 - t) * v0 + t * v1;
      td = (input_p - tfp).norm();
      tr = (1.0 - t) * v[0].radius + t * v[1].radius;
      double tmp = abs(td - tr);
      // logger().debug("edges tmp: {}", tmp);

      if (tmp < min_dist) min_dist = tmp;
    }
  }

  // logger().debug("min_dist after edges loop: {}", min_dist);

  // loop MAT facets
  set<int> near_faces = vertices[min_index].second->faces_;
  for (set<int>::iterator si = near_faces.begin(); si != near_faces.end();
       si++) {
    if (!faces[*si].first) continue;
    NonManifoldMesh_Face mat_f = *faces[*si].second;
    if (mat_f.is_valid_st == false || mat_f.st[0].normal == Vector3::Zero() ||
        mat_f.st[1].normal == Vector3::Zero())
      continue;

    Vector3 v[3], tfp;
    double td;
    for (int i = 0; i < 2; i++) {
      v[0] = mat_f.st[i].v[0];
      v[1] = mat_f.st[i].v[1];
      v[2] = mat_f.st[i].v[2];
      project_point_on_plane(v[0], get_normal(v[0], v[1], v[2]), input_p, tfp,
                             td);
      if (td < min_dist) min_dist = td;
    }
  }

  // logger().debug("min_dist after facets loop: {}", min_dist);
  return min_dist;
}

bool NonManifoldMesh::ValidVertex(int vid) {
  if (vid > vertices.size()) return false;
  return vertices[vid].first;
}

bool NonManifoldMesh::Edge(int vid0, int vid1, int& eid) {
  if (!ValidVertex(vid0) || !ValidVertex(vid1)) return false;

  for (std::set<int>::iterator si = (*vertices[vid0].second).edges_.begin();
       si != (*vertices[vid0].second).edges_.end(); si++) {
    if (edges[*si].first) {
      if (edges[*si].second->HasVertex(vid1)) {
        eid = *si;
        return true;
      }
    }
  }
  return false;
}

bool NonManifoldMesh::Face(const std::set<int>& vset, int& fid) {
  if (vset.size() <= 0) return false;

  for (std::set<int>::iterator si = vset.begin(); si != vset.end(); si++)
    if (!ValidVertex(*si)) return false;

  int vid0 = *(vset.begin());
  for (std::set<int>::iterator si = vertices[vid0].second->faces_.begin();
       si != vertices[vid0].second->faces_.end(); si++)
    if (faces[*si].first) {
      if (faces[*si].second->vertices_ == vset) {
        fid = *si;
        return true;
      }
    }

  return false;
}

void NonManifoldMesh::UpdateCentroid(int fid) {
  if (!faces[fid].first) return;
  if (faces[fid].second->vertices_.size() < 3) {
    log_and_throw(
        "ERROR: we dunno know how to calculate centroid when MAT face has < 3 "
        "vertices");
  }

  faces[fid].second->centroid = Vector3::Zero();
  int count = 0;
  for (const auto& v : faces[fid].second->vertices_) {
    faces[fid].second->centroid += vertices[v].second->pos;
    count++;
  }
  faces[fid].second->centroid /= count;
}

void NonManifoldMesh::UpdateNormal(int fid) {
  if (!faces[fid].first) return;
  if (faces[fid].second->vertices_.size() < 3) {
    log_and_throw(
        "ERROR: we dunno know how to calculate normal when MAT face has < 3 "
        "vertices");
  }

  std::vector<Vector3> vec;
  for (const auto& v : faces[fid].second->vertices_) {
    vec.push_back(vertices[v].second->pos);
  }

  faces[fid].second->normal = (vec[1] - vec[0]).cross(vec[2] - vec[0]);
  faces[fid].second->normal.normalize();
}

void NonManifoldMesh::ComputeFacesNormal() {
  for (int i = 0; i < faces.size(); i++)
    if (faces[i].first) UpdateNormal(i);
}

void NonManifoldMesh::ComputeFacesCentroid() {
  for (int i = 0; i < faces.size(); i++)
    if (faces[i].first) UpdateCentroid(i);
}

void NonManifoldMesh::compute_edges_cone() {
  for (int i = 0; i < edges.size(); i++)
    if (edges[i].first) compute_edge_cone(i);
}

void NonManifoldMesh::compute_edge_cone(int eid) {
  if (!edges[eid].first) return;

  // test validation
  Vector3 c0 = vertices[edges[eid].second->vertices_.first].second->pos;
  Vector3 c1 = vertices[edges[eid].second->vertices_.second].second->pos;
  double r0 = vertices[edges[eid].second->vertices_.first].second->radius;
  double r1 = vertices[edges[eid].second->vertices_.second].second->radius;
  Vector3 c0c1 = c1 - c0;
  double templeng = c0c1.norm() - std::abs(r1 - r0);

  Cone newc(c0, r0, c1, r1);
  edges[eid].second->cone = newc;
  if (newc.type == Cone::Type::INVALID)
    edges[eid].second->is_valid_cone = false;
  else
    edges[eid].second->is_valid_cone = true;
}

void NonManifoldMesh::compute_faces_simple_triangles() {
  for (int i = 0; i < faces.size(); i++)
    if (faces[i].first) compute_face_simple_triangles(i);
}

void NonManifoldMesh::compute_face_simple_triangles(int fid) {
  if (!faces[fid].first) return;

  int k = faces[fid].second->vertices_.size();
  if (k != 3) {
    logger().debug("MAT f: {} has {} vertices", fid, k);
    log_and_throw("ERROR: sorry we don't know how to handle non-tiangle MAT");
  }

  SimpleTriangle st0, st1;
  Vector3 pos[k];
  double radius[k];
  int count = 0;

  for (std::set<int>::iterator si = faces[fid].second->vertices_.begin();
       si != faces[fid].second->vertices_.end(); si++, count++) {
    pos[count] = vertices[*si].second->pos;
    radius[count] = vertices[*si].second->radius;
  }

  int type = get_triangles_from_three_spheres(
      pos[0], radius[0], pos[1], radius[1], pos[2], radius[2], st0, st1);
  faces[fid].second->type = static_cast<NonManifoldMesh_Face::Type>(type);
  if (type == NonManifoldMesh_Face::Type::INVALID) {
    faces[fid].second->is_valid_st = false;
  } else {
    faces[fid].second->st[0] = st0;
    faces[fid].second->st[1] = st1;
    faces[fid].second->is_valid_st = true;
  }
}

void NonManifoldMesh::check_and_store_unthin_tets_in_mat(
    bool is_force_update_neighs) {
  tets_in_mat.clear();
  clear_all_tets();

  // update mat_v_neighbors map
  if (mat_v_neighbors.empty() || is_force_update_neighs) {
    mat_v_neighbors.clear();
    for (const auto face : faces) {
      std::vector<int> f_vs;
      for (const auto v : face.second->vertices_) {
        f_vs.push_back(v);
      }
      if (f_vs.size() != 3) {
        logger().error("mat face has < 3 vertices: {}", f_vs);
        log_and_throw("ERROR");
      }
      for (int i = 0; i < f_vs.size(); i++) {
        mat_v_neighbors[f_vs[i]].insert(f_vs[(i + 1) % 3]);
        mat_v_neighbors[f_vs[i]].insert(f_vs[(i + 2) % 3]);
      }
    }
  }

  // now we have mat_v_neighbors map,
  // for every face with 3 vs, we check if another mat vertex
  // is adjacent to all 3 vs. If so, then tet exits in mat, not thin
  for (const auto face : faces) {
    std::vector<int> f_vs, one_unthin;
    for (const auto v : face.second->vertices_) {
      f_vs.push_back(v);
      one_unthin.push_back(v);
    }
    std::vector<int> intersections;
    set_intersection<int>(mat_v_neighbors.at(f_vs[0]),
                          mat_v_neighbors.at(f_vs[1]),
                          mat_v_neighbors.at(f_vs[2]), intersections);
    std::copy(intersections.begin(), intersections.end(),
              std::back_inserter(one_unthin));

    // if (!intersections.empty()) {
    // 	logger().error("following mat vertices creates a tet: {}",
    // intersections);
    // }
    // std::sort(intersections.begin(), intersections.end());
    // tets_in_mat.insert(intersections);

    for (const auto& v_inter : intersections) {
      bool is_tet_all_face_exit = true;
      std::vector<int> tet = f_vs;
      tet.push_back(v_inter);
      // logger().debug("checking tet: {}", tet);
      // check if every mat face exits
      for (int i = 0; i < 4; i++) {
        for (int j = i + 1; j < 4; j++) {
          for (int k = j + 1; k < 4; k++) {
            std::set<int> vvid;
            vvid.insert(tet[i]);
            vvid.insert(tet[j]);
            vvid.insert(tet[k]);
            int fid;
            if (Face(vvid, fid))  // check if mat face exist
              continue;
            else {
              is_tet_all_face_exit = false;
              // logger().debug("face {} not exist", vvid);
              break;
            }
          }  // k
          if (!is_tet_all_face_exit) break;
        }  // j
        if (!is_tet_all_face_exit) break;
      }  // i

      // save tet if all 4 faces are mat face
      if (is_tet_all_face_exit) {
        std::sort(tet.begin(), tet.end());
        tets_in_mat.insert(tet);
      }
    }  // for intersections
  }

  // creating tets
  for (const auto& tet : tets_in_mat) {
    create_tet(tet);
  }

  logger().debug("[UNTHIN CHECK] tets_in_mat size: {}, tets in MAT {}",
                 tets_in_mat.size(), numTets);
}

void NonManifoldMesh::trace_unthin_mat_faces() {
  // std::queue<int> next_fs_to_check;
  // std::set<int> vs_visited; // mat vertices
  // std::set<int> es_visited; // mat edges
  // std::set<int> fs_visited; // mat faces
  // int prev_euler = 1; // ideal
  // int num_undeleted_fs = 0; // < fs_visited.size()

  if (unthin_trace.empty() && faces.size() > 0)
    next_fs_to_check.push(0);
  else {
    // add last face to propogate
    // next_fs_to_check.push(unthin_trace.back());
  }

  std::vector<int> new_fs_visited;
  while (!next_fs_to_check.empty()) {
    int cur_f = next_fs_to_check.front();
    next_fs_to_check.pop();

    if (fs_visited.find(cur_f) != fs_visited.end()) continue;
    fs_visited.insert(cur_f);

    auto& face = *(faces[cur_f].second);
    if (face.is_deleted) continue;
    unthin_trace.push_back(cur_f);

    // mat face is not deleted
    num_undeleted_fs++;
    new_fs_visited.push_back(cur_f);
    vs_visited.insert(face.vertices_.begin(), face.vertices_.end());
    es_visited.insert(face.edges_.begin(), face.edges_.end());

    // propogate
    for (const auto& e_id : face.edges_) {
      auto& neigh_fs = edges[e_id].second->faces_;
      for (const auto f_id : neigh_fs) {
        next_fs_to_check.push(f_id);
      }
    }

    // check euler
    int euler = vs_visited.size() - es_visited.size() + num_undeleted_fs;
    if (euler != prev_euler) {
      // if (euler > prev_euler) {
      logger().debug("face {} makes euler change {} -> {}", cur_f, prev_euler,
                     euler);
      logger().debug("new_fs_visited: {}", new_fs_visited);
      // found
      unthin_faces.insert(cur_f);
      prev_euler = euler;
      break;
    }
  }

  logger().debug("[Unthin] found {} problematic faces", unthin_faces.size());
}

bool NonManifoldMesh::is_feature_edge_exists(const int v1, const int v2) const {
  if (mat_intf_edges.empty()) return false;

  std::array<int, 2> edge = {{v1, v2}};
  std::sort(edge.begin(), edge.end());
  if (mat_intf_edges.find(edge) != mat_intf_edges.end())
    return true;
  else
    return false;
}

void NonManifoldMesh::insert_new_feature_edge(const int v1, const int v2) {
  std::array<int, 2> edge = {{v1, v2}};
  std::sort(edge.begin(), edge.end());
  mat_intf_edges.insert(edge);

  mat_intf_adjs[v1].insert(v2);
  mat_intf_adjs[v2].insert(v1);
}

// deprecating
void NonManifoldMesh::reload_intf_kd_tree() {
  // mat_intf_adjs.clear();
  // // save mat_intf_edges in another way
  // for (auto const& intf_pair: mat_intf_edges) {
  // 	mat_intf_adjs[intf_pair[0]].insert(intf_pair[1]);
  // 	mat_intf_adjs[intf_pair[1]].insert(intf_pair[0]);
  // }
  intf_kd_tree.reset();
  if (intf_kd_tree == nullptr)
    intf_kd_tree = GEO::NearestNeighborSearch::create(3);
  intf_kd_tree_idx_to_intf_adjs
      .clear();  // idx from intf_kd_points to mat_intf_edges key
  intf_kd_points.clear();
  for (const auto& intf_pair : mat_intf_adjs) {
    int intf_p_idx = intf_pair.first;  // valid_idx
    int intf_p_kd_idx = intf_kd_points.size() / 3;
    for (int i = 0; i < 3; i++) {
      intf_kd_points.push_back(vertices[intf_p_idx].second->pos[i]);
    }
    intf_kd_tree_idx_to_intf_adjs[intf_p_kd_idx] = intf_p_idx;
  }
  intf_kd_tree->set_points(intf_kd_points.size() / 3, intf_kd_points.data());
  logger().error("------ intf_kd_points.size()/3: {}",
                 intf_kd_points.size() / 3);
  logger().error("------- intf_kd_tree->nb_points(): {}",
                 intf_kd_tree->nb_points());
}

///////////////////////////////////////////////////////////////////////////////////////////////
/* # .ma format
 * # all flag_delete might not exist
 * numVertices numEdges numFaces
 * v x y z r flag_type flag_delete
 * e v1 v2 flag_type flag_delete
 * f v1 v2 v3 flag_delete
 */
void MatIO::load_nmm(const std::string& path, NonManifoldMesh& mat) {
  mat.mat_name = getFileName(path, false);
  logger().debug("start loading .ma file: {}", mat.mat_name);
  std::ifstream mastream(path.c_str());
  mat.clear();
  int nv, ne, nf;
  mastream >> nv >> ne >> nf;

  // load MAT vertices
  for (int i = 0; i < nv; i++) {
    char ch;
    double x, y, z, r;
    int type, flag_delete;
    mastream >> ch >> x >> y >> z >> r >> type >> flag_delete;

    Bool_VertexPointer bvp;
    bvp.first = true;
    bvp.second = new NonManifoldMesh_Vertex;
    (*bvp.second).tag = i;
    (*bvp.second).pos[0] = x;
    (*bvp.second).pos[1] = y;
    (*bvp.second).pos[2] = z;
    (*bvp.second).radius = r;
    (*bvp.second).type = static_cast<SphereType>(type);
    (*bvp.second).is_deleted = flag_delete;

    mat.vertices.push_back(bvp);
    mat.numVertices++;
    if (!flag_delete) mat.numVertices_active++;
    mat.max_radius = std::max(mat.max_radius, (*bvp.second).radius);
  }

  // load MAT edges
  for (int i = 0; i < ne; i++) {
    char ch;
    int v1, v2;
    int flag_type;  // is medial feature
                    // (1 - internal, 2 - external, 0 - nothing)
    mastream >> ch;
    mastream >> v1;
    mastream >> v2;
    mastream >> flag_type;
    Bool_EdgePointer bep;
    bep.first = true;
    bep.second = new NonManifoldMesh_Edge;
    (*bep.second).tag = i;
    (*bep.second).vertices_.first = v1;
    (*bep.second).vertices_.second = v2;
    (*mat.vertices[(*bep.second).vertices_.first].second)
        .edges_.insert(mat.edges.size());
    (*mat.vertices[(*bep.second).vertices_.second].second)
        .edges_.insert(mat.edges.size());
    mat.edges.push_back(bep);
    mat.numEdges++;
    mat.numEdges_active++;

    // save medial feature edges
    if (flag_type == 1) {
      std::array<int, 2> edge = {{v1, v2}};
      std::sort(edge.begin(), edge.end());
      mat.mat_intf_edges.insert(edge);
    } else if (flag_type == 2) {
      std::array<int, 2> edge = {{v1, v2}};
      std::sort(edge.begin(), edge.end());
      mat.mat_extf_edges.insert(edge);
    }
  }

  // load MAT faces
  for (int i = 0; i < nf; i++) {
    char ch;
    int vid[3];
    int eid[3];
    // int is_deleted;
    // mastream >> ch >> vid[0] >> vid[1] >> vid[2] >> is_deleted;
    mastream >> ch >> vid[0] >> vid[1] >> vid[2];

    logger().debug("mat f {} has vertices ({},{},{})", vid[0], vid[1], vid[2]);
    if (vid[0] == vid[1] || vid[0] == vid[2] || vid[1] == vid[2]) {
      logger().error("MAT face {} is degenerated, skip", i);
      continue;
    }

    Bool_FacePointer bfp;
    bfp.first = true;
    bfp.second = new NonManifoldMesh_Face;
    (*bfp.second).tag = i;
    for (int v = 0; v < 3; v++) {
      (*bfp.second).vertices_.insert(vid[v]);
    }
    // (*bfp.second).is_deleted = is_deleted == 1 ? true : false;

    if (mat.Edge(vid[0], vid[1], eid[0])) (*bfp.second).edges_.insert(eid[0]);
    if (mat.Edge(vid[0], vid[2], eid[1])) (*bfp.second).edges_.insert(eid[1]);
    if (mat.Edge(vid[1], vid[2], eid[2])) (*bfp.second).edges_.insert(eid[2]);
    mat.vertices[vid[0]].second->faces_.insert(mat.faces.size());
    mat.vertices[vid[1]].second->faces_.insert(mat.faces.size());
    mat.vertices[vid[2]].second->faces_.insert(mat.faces.size());
    mat.edges[eid[0]].second->faces_.insert(mat.faces.size());
    mat.edges[eid[1]].second->faces_.insert(mat.faces.size());
    mat.edges[eid[2]].second->faces_.insert(mat.faces.size());
    mat.faces.push_back(bfp);
    mat.numFaces++;
    mat.numFaces_active++;
  }

  mat.computebb();
  mat.ComputeFacesCentroid();
  mat.ComputeFacesNormal();
  mat.compute_edges_cone();
  mat.compute_faces_simple_triangles();
  mat.print_info();
}

// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
//
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void MatIO::export_ma_given(const std::string& maname,
                            const std::vector<Vector4>& mat_vertices,
                            const std::vector<std::array<int, 2>>& mat_edges,
                            const std::vector<std::array<int, 3>>& mat_faces,
                            bool is_use_given_name) {
  std::string ma_name_full = maname;
  if (!is_use_given_name)
    ma_name_full = "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ma";

  std::ofstream fout;
  fout.open(ma_name_full, std::ofstream::out | std::ofstream::app);  //   append
  fout << mat_vertices.size() << " " << mat_edges.size() << " "
       << mat_faces.size() << std::endl;

  // save vertices
  for (int i = 0; i < mat_vertices.size(); i++) {
    const auto& mat_v = mat_vertices.at(i);
    fout << "v " << std::setiosflags(std::ios::fixed) << std::setprecision(15)
         << mat_v[0] << " " << mat_v[1] << " " << mat_v[2] << " " << mat_v[3];
    fout << std::endl;
  }

  //  save edges
  for (int i = 0; i < mat_edges.size(); i++) {
    const auto& mat_e = mat_edges[i];
    fout << "e " << mat_e[0] << " " << mat_e[1];
    fout << std::endl;
  }

  // save faces
  for (int i = 0; i < mat_faces.size(); i++) {
    const auto& mat_f = mat_faces[i];
    fout << "f";
    for (uint v = 0; v < 3; v++) fout << " " << mat_f[v];
    fout << std::endl;
  }
  fout.close();

  printf("saved mat at: %s \n", ma_name_full.c_str());
}

// helper function for export_ma_clean() and export_ma_ply()
void MatIO::get_mat_clean(const NonManifoldMesh& mat,
                          std::vector<Vector4>& vertices,
                          std::vector<std::array<int, 2>>& edges,
                          std::vector<std::array<int, 3>>& faces) {
  std::map<int, int> map_vertices;  // mat vertex tag to new
  vertices.clear();
  edges.clear();
  faces.clear();

  auto get_vertex_mapped_id = [&](const int vid) {
    if (map_vertices.find(vid) == map_vertices.end()) {
      // add a new vertex
      map_vertices[vid] = map_vertices.size();
    }
    return map_vertices.at(vid);
  };

  // faces
  for (int f = 0; f < mat.faces.size(); f++) {
    const auto& face = *mat.faces[f].second;
    if (face.is_deleted) continue;
    std::array<int, 3> one_f;
    for (uint j = 0; j < 3; j++) {
      int vid = get_vertex_mapped_id(*std::next(face.vertices_.begin(), j));
      one_f[j] = vid;
    }
    faces.push_back(one_f);
  }
  printf("faces: %d \n", faces.size());

  // edges
  for (int e = 0; e < mat.edges.size(); e++) {
    const auto& edge = *mat.edges[e].second;
    if (edge.is_deleted) continue;
    // if (edge.faces_.empty()) continue;
    int vid1 = get_vertex_mapped_id(edge.vertices_.first);
    int vid2 = get_vertex_mapped_id(edge.vertices_.second);
    edges.push_back({{vid1, vid2}});
  }
  printf("edges: %d \n", edges.size());

  // vertices
  // save from map_vertices, to avoid (0,0,0,0) in .ma file
  vertices.resize(map_vertices.size());
  for (const auto& v_pair : map_vertices) {
    int old_vid = v_pair.first;
    int new_vid = v_pair.second;
    const auto& mat_v = *mat.vertices[old_vid].second;
    vertices[new_vid] =
        Vector4(mat_v.pos[0], mat_v.pos[1], mat_v.pos[2], mat_v.radius);
  }

  printf("vertcies: %d \n", vertices.size());
}

// remove deleted edges/faces
// should save the same result as export_ma_ply but with different format
//
// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
/* # .ma format
 * numVertices numEdges numFaces
 * v x y z r
 * e v1 v2
 * f v1 v2 v3
 */
void MatIO::export_ma_clean(const std::string& maname,
                            const NonManifoldMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ma";
  printf("start saving mat .m file: %s \n", ma_name_full.c_str());

  std::vector<Vector4> vertices;
  std::vector<std::array<int, 2>> edges;
  std::vector<std::array<int, 3>> faces;
  get_mat_clean(mat, vertices, edges, faces);

  // export
  export_ma_given(maname, vertices, edges, faces);
}

// this export is matching blender addon
// https://github.com/songshibo/blender-mat-addon
/* # .ma format
 * # all flag_delete might not exist
 * numVertices numEdges numFaces
 * v x y z r flag_type flag_delete
 * e v1 v2 flag_type flag_delete
 * f v1 v2 v3 flag_delete
 */
void MatIO::export_nmm(const std::string& maname, const NonManifoldMesh& mat) {
  // std::string ma_name_full =
  //     "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ma";
  std::string ma_name_full = "../out/mat/mat_" + maname + ".ma";
  logger().debug("start saving .ma file: {}", ma_name_full);
  const auto& mat_intf_edges = mat.mat_intf_edges;
  const auto& mat_extf_edges = mat.mat_extf_edges;

  std::ofstream fout(ma_name_full);
  fout << mat.numVertices << " " << mat.numEdges_active << " "
       << mat.numFaces_active << std::endl;

  // we might have redundant vertices
  // but we need those indices, so store it as delete flag = true
  for (int i = 0; i < mat.vertices.size(); i++) {
    Vector3 pos = mat.vertices[i].second->pos;
    fout << "v " << setiosflags(ios::fixed) << setprecision(15) << pos[0] << " "
         << pos[1] << " " << pos[2] << " " << mat.vertices[i].second->radius
         << " " << int(mat.vertices[i].second->type);

    // store this for not showing redundant vertices
    if (mat.vertices[i].second->is_deleted) {
      fout << " 1";
    } else {
      fout << " 0";
    }
    fout << std::endl;
  }

  // do not save edge if deleted
  for (int i = 0; i < mat.edges.size(); i++) {
    if (mat.edges[i].second->is_deleted) continue;
    std::array<int, 2> edge = {{mat.edges[i].second->vertices_.first,
                                mat.edges[i].second->vertices_.second}};
    std::sort(edge.begin(), edge.end());
    fout << "e " << edge[0] << " " << edge[1];
    // is this edge a feature edge
    if (mat_intf_edges.find(edge) != mat_intf_edges.end()) {
      fout << " " << 1;  // internal
    } else if (mat_extf_edges.find(edge) != mat_extf_edges.end()) {
      fout << " " << 2;  // external
    } else {
      fout << " " << 0;
    }
    fout << std::endl;
  }

  // do not save face if deleted
  for (int i = 0; i < mat.faces.size(); i++) {
    if (mat.faces[i].second->is_deleted) continue;
    fout << "f";
    for (auto si = mat.faces[i].second->vertices_.cbegin();
         si != mat.faces[i].second->vertices_.cend(); si++)
      fout << " " << *si;

    fout << std::endl;
  }
  fout.close();
}

// some vertices/faces are deleted
void MatIO::write_nmm_ply(const std::string& maname,
                          const NonManifoldMesh& mat) {
  // std::string ma_name_full =
  // "../out/mat/mat_" + maname + "_" + get_timestamp() + ".ply";
  std::string ma_name_full = "../out/mat/mat_" + maname + ".ply";
  logger().debug("start saving mat .ply file: {}", ma_name_full);

  std::map<int, int> map_vertices;  // mat vertex tag to new ply
  int num_vertices = 0;
  int num_edges = 0;
  int num_faces = 0;
  std::vector<Vector3> vertices;
  std::vector<std::array<int, 2>> edges;
  std::vector<std::array<int, 3>> faces;

  auto get_vertex_mapped_id = [&](const int vid) {
    if (map_vertices.find(vid) == map_vertices.end()) {
      map_vertices[vid] = num_vertices;
      num_vertices++;
    }
    return map_vertices.at(vid);
  };

  for (int f = 0; f < mat.faces.size(); f++) {
    const auto& face = *(mat.faces[f].second);
    if (face.is_deleted) continue;
    num_faces++;
    int j = 0;
    std::array<int, 3> one_f;
    for (auto si = face.vertices_.cbegin(); si != face.vertices_.cend();
         si++, j++) {
      int vid = get_vertex_mapped_id(*si);
      one_f[j] = vid;
    }
    faces.push_back(one_f);
  }
  logger().debug("faces: {}", num_faces);

  for (int e = 0; e < mat.edges.size(); e++) {
    const auto& edge = *(mat.edges[e].second);
    if (edge.is_deleted) continue;
    num_edges++;
    int vid1 = get_vertex_mapped_id(edge.vertices_.first);
    int vid2 = get_vertex_mapped_id(edge.vertices_.second);
    edges.push_back({{vid1, vid2}});
  }
  logger().debug("edges: {}", num_edges);

  vertices.resize(map_vertices.size());
  for (int v = 0; v < mat.vertices.size(); v++) {
    auto& vertex = *(mat.vertices[v].second);
    if (vertex.is_deleted) continue;
    int vid = get_vertex_mapped_id(vertex.tag);
    if (vid == vertices.size()) {
      vertices.push_back(vertex.pos);
    } else {
      vertices[vid] = vertex.pos;
    }
  }
  logger().debug("vertcies: {}", num_vertices);

  MatrixXs V(num_vertices, 3);
  Eigen::MatrixXi E(num_edges, 2);
  Eigen::MatrixXi F(num_faces, 3);
  for (int v = 0; v < num_vertices; v++) {
    V(v, 0) = vertices[v][0];
    V(v, 1) = vertices[v][1];
    V(v, 2) = vertices[v][2];
  }
  for (int e = 0; e < num_edges; e++) {
    E(e, 0) = edges[e][0];
    E(e, 1) = edges[e][1];
  }
  for (int f = 0; f < num_faces; f++) {
    F(f, 0) = faces[f][0];
    F(f, 1) = faces[f][1];
    F(f, 2) = faces[f][2];
  }
  igl::writePLY(ma_name_full, V, F, E);
}

void MatIO::export_nmm_vs(const std::string& maname,
                          const NonManifoldMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_" + maname + "_" + get_timestamp() + ".r";
  logger().debug("start saving mat vertices to .r file: {}", ma_name_full);

  std::ofstream fout(ma_name_full);
  fout << mat.numVertices << std::endl;

  for (int i = 0; i < mat.vertices.size(); i++) {
    Vector3 pos = mat.vertices[i].second->pos;
    fout << setiosflags(ios::fixed) << setprecision(15) << pos[0] << " "
         << pos[1] << " " << pos[2] << " " << mat.vertices[i].second->radius
         << std::endl;
  }
  fout.close();
}

void MatIO::export_nmm_tets(const std::string& maname,
                            const NonManifoldMesh& mat) {
  std::string ma_name_full =
      "../out/mat/mat_tets_" + maname + "_" + get_timestamp() + ".ma";
  logger().debug("start saving mat tets to .ma file: {}", ma_name_full);

  std::vector<int> tet_vs_vec;
  std::vector<Vector3i> tet_faces;
  for (int t = 0; t < mat.tets.size(); t++) {
    const auto& tet_vs = mat.tets[t].vertices_;
    if (tet_vs.size() != 4) log_and_throw("ERROR: tet_vs size not 4");
    tet_vs_vec.clear();
    for (const auto& vs : tet_vs) tet_vs_vec.push_back(vs);
    tet_faces.push_back(Vector3i(tet_vs_vec[0], tet_vs_vec[1], tet_vs_vec[2]));
    tet_faces.push_back(Vector3i(tet_vs_vec[0], tet_vs_vec[1], tet_vs_vec[3]));
    tet_faces.push_back(Vector3i(tet_vs_vec[0], tet_vs_vec[2], tet_vs_vec[3]));
    tet_faces.push_back(Vector3i(tet_vs_vec[1], tet_vs_vec[2], tet_vs_vec[3]));
  }

  std::ofstream fout(ma_name_full);
  fout << mat.numVertices << " " << 0 << " " << tet_faces.size() << std::endl;

  // we might have redundant vertices
  // but we need those indices, so store it as delete flag = true
  for (int i = 0; i < mat.vertices.size(); i++) {
    Vector3 pos = mat.vertices[i].second->pos;
    fout << "v " << setiosflags(ios::fixed) << setprecision(15) << pos[0] << " "
         << pos[1] << " " << pos[2] << " " << mat.vertices[i].second->radius
         << " " << int(mat.vertices[i].second->type);

    // store this for not showing redundant vertices
    if (mat.vertices[i].second->is_deleted) {
      fout << " 1";
    } else {
      fout << " 0";
    }
    fout << std::endl;
  }

  for (int f = 0; f < tet_faces.size(); f++) {
    fout << "f";
    for (int i = 0; i < 3; i++) {
      fout << " " << tet_faces[f][i];
    };
    fout << std::endl;
  }

  fout.close();
}

}  // namespace matfp