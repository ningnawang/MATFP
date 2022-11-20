// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/Triangulation.h"

#include <geogram/basic/progress.h>
#include <geogram/mesh/triangle_intersection.h>
#include <igl/parallel_for.h>

#include <cassert>
#include <mutex>

#include "matfp/Common.h"
#include "matfp/InscribedSpheres.h"
#include "matfp/LFS.h"
#include "matfp/MeshProcessor.h"
#include "matfp/WindingFilter.h"
#include "matfp/geogram/RPD.h"
#define USE_MVERTEX false
#define USE_MVERTEX_FOR_MAT true

namespace matfp {
//////////////////////////////////////////////////////////////////////
// Helper functions
/////////////////////////////////////////////////////////////////////
//
// will update sp_edges and sp_corners
void store_external_feature_spheres(const std::vector<double> &seed_points,
                                    std::vector<SpecialEdge> &sp_edges,
                                    std::vector<SpecialCorner> &sp_corners,
                                    std::vector<MVertex> &all_medial_spheres,
                                    std::set<std::array<int, 2>> &se_spheres) {
  se_spheres.clear();
  std::vector<Vector3> se_normal;
  std::vector<int> se_ref_fs_ids;        /* matching GEO::Mesh*/
  std::map<int, int> map_seed_2_mat;     // seed_idx -> mat.tag
  std::map<int, int> corners_seed_2_sp;  // seed_idx -> SpecialCorner::id
  for (const auto &sp_c : sp_corners) {
    corners_seed_2_sp[sp_c.seed_id] = sp_c.id;
  }

  logger().debug("start storing external feature spheres...");
  for (auto &one_sp_edge : sp_edges) {
    if (one_sp_edge.type == EdgeType::CCE) continue;
    const auto &one_se = one_sp_edge.vs_seed_ids;
    const auto &se_ns = one_sp_edge.adj_ref_normals;
    const auto &se_ref_fs = one_sp_edge.adj_ref_fs_pair;
    se_normal.clear();
    std::copy(se_ns.begin(), se_ns.end(), std::back_inserter(se_normal));
    se_ref_fs_ids.clear();
    std::copy(se_ref_fs.begin(), se_ref_fs.end(),
              std::back_inserter(se_ref_fs_ids));
    for (int i = 0; i < 2; i++) {
      const int seed_id = one_se[i];
      // create new MVertex
      if (map_seed_2_mat.find(seed_id) == map_seed_2_mat.end()) {
        MVertex f_mat_p(
            all_medial_spheres.size(), seed_points[seed_id * 3],
            seed_points[seed_id * 3 + 1], seed_points[seed_id * 3 + 2],
            SCALAR_FEATURE_SQ_RADIUS /*sq_radius*/, SphereType::T_1_N);
        f_mat_p.is_on_s_edge = true;
        f_mat_p.dilate_sphere_radius();
        // update tangent planes
        f_mat_p.add_tan_planes_for_feature_sphere(se_normal, se_ref_fs_ids);
        all_medial_spheres.push_back(f_mat_p);
        map_seed_2_mat[seed_id] = f_mat_p.tag;
      }
      if (corners_seed_2_sp.find(seed_id) != corners_seed_2_sp.end()) {
        // check corners
        // save tangent planes for corners
        const int c_tag = map_seed_2_mat.at(seed_id);
        MVertex &c_mat_p = all_medial_spheres.at(c_tag);
        c_mat_p.add_tan_planes_for_feature_sphere(se_normal, se_ref_fs_ids);
        c_mat_p.is_on_corner = true;
        // update special corners
        int sp_cid = corners_seed_2_sp.at(seed_id);
        sp_corners[sp_cid].mat_tag = c_tag;
      }
    }  // for i < 2

    // add sharp edge neighboring info to all_medial_spheres
    // for later detecting invalid MAT edges
    const int s1_tag = map_seed_2_mat.at(one_se[0]);
    const int s2_tag = map_seed_2_mat.at(one_se[1]);
    all_medial_spheres.at(s1_tag).se_adj_se.insert(s2_tag);
    all_medial_spheres.at(s2_tag).se_adj_se.insert(s1_tag);

    std::array<int, 2> one_pair = {{s1_tag, s2_tag}};
    std::sort(one_pair.begin(), one_pair.end());
    se_spheres.insert(one_pair);
    one_sp_edge.mat_se_tags = one_pair;
  }  // for sp_edges

  // corner sphere: sanity check
  for (const auto &c_pair : corners_seed_2_sp) {
    int c_seed_id = c_pair.first;
    const int c_tag = map_seed_2_mat.at(c_seed_id);
    auto &c_mat_p = all_medial_spheres.at(c_tag);
    if (c_mat_p.tan_planes.size() <= 2) {
      logger().error(
          "corners: mat_p {} has tan_planes {} <= 2, please double check if "
          "this is a real corner!",
          c_mat_p.tag, c_mat_p.tan_planes.size());
      c_mat_p.print_info();
      log_and_throw("ERROR");
    }
  }

  logger().debug("added {} feature mat points, {} se_spheres pairs",
                 map_seed_2_mat.size(), se_spheres.size());
}

//////////////////////////////////////////////////////////////////////
// Main functions
/////////////////////////////////////////////////////////////////////

void update_dt(ThreeDimensionalShape &shape3D, bool is_debug) {
  if (is_debug) logger().debug("calling update_dt ...");
  igl::Timer timer;
  const std::vector<double> &seed_points = shape3D.seed_points;
  const std::vector<Vector3> &seed_normals = shape3D.seed_normals_v2;
  const std::set<int> &feature_points = shape3D.feature_points;
  const std::set<int> &seed_is_deleted = shape3D.seed_is_deleted;

  // bbox
  bool &is_bbox_stored = shape3D.params.is_bbox_stored;
  std::vector<Vector3> &bb_points = shape3D.params.bb_points;

  // will be updated
  std::vector<MVertex> &all_medial_spheres = shape3D.all_medial_spheres;
  std::vector<SpecialCorner> &sp_corners = shape3D.sp_corners;
  std::vector<SpecialEdge> &sp_edges = shape3D.sp_edges;
  DelaunayTriangulation &dt = shape3D.dt;

  const int dim = 3;
  std::vector<double> dt_vs, dt_vs_normals;

  // sanity check
  if (seed_points.size() != seed_normals.size() * 3) {
    logger().debug("seed_points size: {}, seed_normals size: {}",
                   seed_points.size() / 3, seed_normals.size());
    log_and_throw("Error");
  }

  // // I: initial Delaunay vertices
  // step 1: all seed points
  // 		   (include both feature and non-feature points)
  // NOTE: we do not handle seed_is_deleted here
  //		 because we want dt vertices has tags that is mapping to
  // seed_points
  int nb_pts = seed_points.size() / dim;
  for (int i = 0; i < nb_pts; i++) {
    for (int j = 0; j < dim; j++) {
      dt_vs.push_back(seed_points[i * 3 + j]);
      dt_vs_normals.push_back(seed_normals[i][j]);
    }
  }
  // step 2: add 8 bbox
  if (!is_bbox_stored) {
    Vector3 min, max;
    get_bb_corners(shape3D.params, seed_points, min, max);
    for (int i = 0; i < 8; i++) {
      std::array<double, 3> p;
      std::bitset<sizeof(int) * 8> flag(i);
      for (int j = 0; j < 3; j++) {
        if (flag.test(j))
          p[j] = max[j];
        else
          p[j] = min[j];
      }
      bb_points.push_back(Vector3(p[0], p[1], p[2]));
    }
    is_bbox_stored = true;
  }
  if (is_debug) logger().info("adding 8 bbox ...");
  for (int i = 0; i < 8; i++) {
    for (int j = 0; j < 3; j++) {
      dt_vs.push_back(bb_points[i][j]);
      dt_vs_normals.push_back(0.);  // add (0,0,0) as normal for 8 bbox
    }
  }
  if (is_debug)
    logger().info("Initial DT points: {} = total {} + bbox 8", dt_vs.size() / 3,
                  seed_points.size() / 3);
  if (dt_vs.size() != dt_vs_normals.size())
    log_and_throw("ERROR: dt_vs.size() != dt_vs_normals.size()");

  // step 3: calculate weight for non-feature mat points
  // would store non-feature mat points
  matfp::generate_DT_CGAL(dt_vs, dt_vs_normals, seed_is_deleted, feature_points,
                          dt, all_medial_spheres, is_debug);

  timer.start();
  // step 4: in/out filter non-feature mat points
  matfp::inout_filter_raw(all_medial_spheres, shape3D.sf_mesh_wrapper.VI,
                          shape3D.sf_mesh_wrapper.FI, shape3D.winding_num,
                          shape3D.is_outside, shape3D.dt);
  if (is_debug)
    logger().info("winding filter took {}s", timer.getElapsedTimeInSec());

  // step 5: collect tangent info of non-feature mat points
  collect_and_update_DT_feature_tets(shape3D.dt, seed_points, seed_normals,
                                     all_medial_spheres);

  // step 6: store feature mat points and the se_spheres pairs
  // will update sp_edges and sp_corners
  store_external_feature_spheres(seed_points, sp_edges, sp_corners,
                                 all_medial_spheres, shape3D.se_spheres);

  // step 7: store bbox spheres to avoid handling infinite RT tets
  // that might impact our MAT construction
  // used in generate_RT_for_dual_PD()
  if (is_debug) logger().info("adding 8 bbox spheres ...");
  for (int i = 0; i < 8; i++) {
    const Vector3 &bb_point = bb_points[i];
    // add new MVertex
    int all_idx = all_medial_spheres.size();
    MVertex mat_bbox(all_idx, bb_point[0], bb_point[1], bb_point[2],
                     SCALAR_FEATURE_SQ_RADIUS);
    mat_bbox.is_on_bbox = true;
    mat_bbox.is_outside = true;
    all_medial_spheres.push_back(mat_bbox);
  }

  // step 8: add offset to inside spheres' radii for avoiding degeneracy
  // adding here instead of generate_RT_CGAL_MVertex()
  // is because function prepare_all_medial_spheres()
  // takes updated radius into account when deleting spheres
  for (auto &mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) {
      continue;
    }
    mat_p.dilate_sphere_radius();
  }
  if (is_debug) logger().debug("done update_dt");
}

// and call remove_medial_spheres_too_close_to_sharp_edges() ahead
void generate_RT_for_dual_PD(std::vector<MVertex> &all_medial_spheres,
                             std::vector<int> &valid_medial_spheres,
                             const std::vector<double> &seed_points,
                             const std::set<int> &feature_points,
                             const std::vector<Vector3> &bb_points,
                             RegularTriangulationNN &rt,
                             bool is_using_dilated_radius, bool is_debug) {
  if (is_debug)
    logger().debug(
        "----- generating Regular Triangulation as dual of Power Diagram");
  rt.clean();

  // add inside medial spheres
  // including sample medial spheres
  for (int all_idx = 0; all_idx < all_medial_spheres.size(); all_idx++) {
    MVertex &mat_p = all_medial_spheres[all_idx];
    if (mat_p.is_deleted) continue;
    // add bbox points to avoid handling infinite RT tets
    // that might impact our MAT construction
    // if is_on_bbox= true then is_outside must be true
    if (mat_p.is_outside && !mat_p.is_on_bbox) continue;
    double sq_radius = mat_p.sq_radius;
    if (is_using_dilated_radius) {
      // update dilated radius before calling
      mat_p.dilate_sphere_radius();
      sq_radius = mat_p.sq_radius_dilated;
    }
    Point_rt p(mat_p.pos[0], mat_p.pos[1], mat_p.pos[2]);
    Weight weight = sq_radius;
    Vertex_handle_rt vh;
    vh = rt.insert(Weighted_point(p, weight));
    if (vh == nullptr) continue;  // THIS IS IMPORTANT
    // update RT
    vh->info().all_tag = mat_p.tag;
    vh->info().pos = mat_p.pos;
    vh->info().radius = std::sqrt(mat_p.sq_radius);
    vh->info().sq_radius = mat_p.sq_radius;
    vh->info().is_feature = mat_p.is_a_feature_sphere();
    vh->info().is_corner = mat_p.is_on_corner;
    vh->info().is_outside = mat_p.is_outside;
    vh->info().is_on_bbox = mat_p.is_on_bbox;
  }

  // clear valid_idx, will assign new after RT
  // NOTE: we add 8 bbox spheres, so
  // valid_medial_spheres.size() + 8 == rt.number_of_vertices()
  valid_medial_spheres.clear();
  for (int all_idx = 0; all_idx < all_medial_spheres.size(); all_idx++) {
    MVertex &mat_p = all_medial_spheres[all_idx];
    mat_p.valid_idx = -1;
  }
  for (Finite_vertices_iterator_rt vit = rt.finite_vertices_begin();
       vit != rt.finite_vertices_end(); vit++) {
    // // is_on_bbox also outside, but we store its valid_idx
    // if (vit->info().is_outside && !vit->info().is_on_bbox) {
    //   logger().error("vit is outside but not on bbox, wrong!!!!");
    //   log_and_throw("ERROR");
    // }
    // do not store valid_idx for sphere on bbox
    if (vit->info().is_outside || vit->info().is_on_bbox) continue;
    // https://stackoverflow.com/a/4206085
    // in CGAL triangulation, iterator = handle
    int valid_idx = valid_medial_spheres.size();
    Vertex_handle_rt vh = vit;
    vh->info().tag = valid_idx;
    rt.set_tag_to_vh(valid_idx, vh);
    int all_idx = vh->info().all_tag;
    // store valid_idx, already set to -1 earlier
    all_medial_spheres[all_idx].valid_idx = valid_idx;
    valid_medial_spheres.push_back(all_idx);
  }
  if (rt.number_of_vertices() != valid_medial_spheres.size() + 8) {
    logger().error(
        "ERROR: valid_medial_spheres.size() + 8: {} != "
        "rt.number_of_vertices(): {}",
        valid_medial_spheres.size() + 8, rt.number_of_vertices());
    log_and_throw("ERROR");
  }
  rt.set_nb_vertices(rt.number_of_vertices());
  rt.print_info();

  // sanity check
  for (int valid_idx = 0; valid_idx < valid_medial_spheres.size();
       valid_idx++) {
    const int all_idx = valid_medial_spheres[valid_idx];
    if (all_idx == -1) {
      logger().debug("checkin valid_idx: {}", valid_idx);
      logger().debug("checkin all_idx: {} cannot be -1!!!", all_idx);
      log_and_throw("ERROR");
    };
  }
}

void get_RT_tet_dual_info(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const AABBWrapper &aabb_wrapper,
                          RegularTriangulationNN &rt) {
  // fetching dual info of RT tets
  logger().debug("[RT dual] fetching dual info of RT tets ...");
  const double evenlop_eps = 0.01;
  std::vector<Vector3> in_out_points;
  std::vector<bool> are_outside;  // store is_outside
  int tid = 0;
  rt.rt_ts_info.clear();
  for (Finite_cells_iterator_rt fci = rt.finite_cells_begin();
       fci != rt.finite_cells_end(); fci++) {
    RTI rti;
    for (int i = 0; i < 4; i++) {
      rti.valid_vs_tags[i] = fci->vertex(i)->info().tag;
      rti.all_vs_tags[i] = fci->vertex(i)->info().all_tag;
      rti.is_contain_corner =
          rti.is_contain_corner || fci->vertex(i)->info().is_corner;
    }
    std::sort(rti.valid_vs_tags.begin(), rti.valid_vs_tags.end());
    std::sort(rti.all_vs_tags.begin(), rti.all_vs_tags.end());
    RegularTriangulationNN::Bare_point bp =
        rt.dual(fci);  // dunno why we do not need *fci here
    rti.center = Vector3(bp[0], bp[1], bp[2]);
    // rti.dist2sf_dual_center =
    //     std::sqrt(aabb_wrapper.get_sq_dist_to_sf(rti.center));
    rti.sq_radius = bp[3];
    // save id only when we save it
    rti.id = tid++;
    rt.rt_ts_info.push_back(rti);
    in_out_points.push_back(rti.center);
    fci->info() = rti;  // save to RT
  }
  // store inside/outside of dual RT tets
  inout_filter_vector_points(in_out_points, V, F, are_outside);
  for (tid = 0; tid < rt.rt_ts_info.size(); tid++) {
    RTI &rti = rt.rt_ts_info[tid];
    rti.is_dual_point_outside = are_outside[tid];
    // // all inside if within narrow band
    // if (rti.dist2sf_dual_center <= evenlop_eps)
    //   rti.is_dual_point_outside = false;
  }
  logger().debug("[RT dual] fetched dual info of RT tets: {}",
                 rt.rt_ts_info.size());
}

void get_RT_face_dual_info(const AABBWrapper &aabb_wrapper,
                           RegularTriangulationNN &rt) {
  // fetching dual info of RT faces
  // a facet is given by a cell and an index
  // (the facet i of a cell c is the facet of c that is opposite to the vertex
  // with index i)
  logger().debug("[RT dual] fetching dual info of RT faces ...");
  int fid = 0;
  for (Finite_facets_iterator_rt fft = rt.finite_facets_begin();
       fft != rt.finite_facets_end(); fft++) {
    RFI rfi;
    // fetch 3 vertices info
    for (int i = 0; i < 3; i++) {
      rfi.valid_vs_tags[i] =
          fft->first->vertex((fft->second + i + 1) % 4)->info().tag;
      rfi.all_vs_tags[i] =
          fft->first->vertex((fft->second + i + 1) % 4)->info().all_tag;
    }
    std::sort(rfi.valid_vs_tags.begin(), rfi.valid_vs_tags.end());
    std::sort(rfi.all_vs_tags.begin(), rfi.all_vs_tags.end());

    // fetch dual segment info
    const auto &cell2 = fft->first->neighbor(int(fft->second));
    rfi.adj_tets = {{fft->first->info().id, cell2->info().id}};

    // check inifinity, do not save
    // we've add 8 bbox sphere during RT, so no face of medial spheres
    // on boundary and adjacent to infinite tet
    if (rt.is_infinite(fft->first) || rt.is_infinite(cell2)) {
      // logger().error("RT face has any tet infinite, rfi.all_vs_tags {}",
      //                rfi.all_vs_tags);
      // rfi.print_info();
      // if (rfi.adj_tets[0] == -1) {
      //   rt.rt_ts_info[rfi.adj_tets[1]].print_info();
      // } else {
      //   rt.rt_ts_info[rfi.adj_tets[0]].print_info();
      // }
      continue;
    }
    // both adjacent tets are finite
    // logger().debug("rfi.adj_tets: {}", rfi.adj_tets);
    rfi.dual_segment = {{rt.rt_ts_info[rfi.adj_tets[0]].center,
                         rt.rt_ts_info[rfi.adj_tets[1]].center}};
    // calcualte the intersection of dual_segment and surface
    // compute the intersection point nearest to the extremity
    // rfi.dual_segment[0]
    rfi.is_dual_seg_intersect = aabb_wrapper.sf_segment_nearest_intersection(
        rfi.dual_segment[0], rfi.dual_segment[1], rfi.dual_intersections[0]);
    if (rfi.is_dual_seg_intersect) {
      // compute the intersection point nearest to the extremity
      // rfi.dual_segment[1] NOTE: changing extremity may give us different
      // boolean result
      bool is_intersect_2 = aabb_wrapper.sf_segment_nearest_intersection(
          rfi.dual_segment[1], rfi.dual_segment[0], rfi.dual_intersections[1]);
      rfi.is_dual_seg_intersect = rfi.is_dual_seg_intersect && is_intersect_2;
      // logger().debug(
      //     "rfi {} dual seg has two intersections ({},{},{}), ({},{},{})",
      //     rfi.id, rfi.dual_intersections[0][0],
      //     rfi.dual_intersections[0][1], rfi.dual_intersections[0][2],
      //     rfi.dual_intersections[1][0], rfi.dual_intersections[1][1],
      //     rfi.dual_intersections[1][2]);
    }
    if (!rfi.is_dual_seg_intersect) {
      rfi.dual_intersections = {{Vector3::Zero(), Vector3::Zero()}};
    }

    // save id only when we save it
    rfi.id = fid++;
    rt.rt_fs_info.push_back(rfi);

    // save RTI::adj_facets
    for (const auto &tid : rfi.adj_tets) {
      auto &rti = rt.rt_ts_info[tid];
      rti.adj_facets.insert(rfi.id);
    }
  }  // for Finite_facets_iterator_rt

  // update rfi.dist_dual_intersections
  for (fid = 0; fid < rt.rt_fs_info.size(); fid++) {
    RFI &rfi = rt.rt_fs_info[fid];
    rfi.is_dual_seg_vs_outside = {
        {rt.rt_ts_info[rfi.adj_tets[0]].is_dual_point_outside,
         rt.rt_ts_info[rfi.adj_tets[1]].is_dual_point_outside}};
    // update the distance of intersection segment
    // (rfi.dist_dual_intersections) between dual_segement and the surface
    if (!rfi.is_dual_seg_vs_outside[0] && !rfi.is_dual_seg_vs_outside[1]) {
      // if both inside
      rfi.dual_intersections = rfi.dual_segment;
    } else if (!rfi.is_dual_seg_vs_outside[0] &&
               rfi.is_dual_seg_vs_outside[1]) {
      // if 1st insdie, 2nd outside
      rfi.dual_intersections = {
          {rfi.dual_segment[0], rfi.dual_intersections[1]}};
    } else if (rfi.is_dual_seg_vs_outside[0] &&
               !rfi.is_dual_seg_vs_outside[1]) {
      // if 1st outside, 2nd inside
      rfi.dual_intersections = {
          {rfi.dual_intersections[0], rfi.dual_segment[1]}};
    } else {
      // if both outside, no change
      // so that dual_seg distance won't be 0
    }
    rfi.dist_dual_intersections = get_distance_between_two_vectors(
        rfi.dual_intersections[0], rfi.dual_intersections[1]);
  }  // for rt.rt_fs_info
  logger().debug("[RT dual] fetched dual info of RT faces: {}",
                 rt.rt_fs_info.size());
}

// helper function for get_RT_edge_dual_info
void get_RT_edge_dual_poly_sf_intersections(const GEO::Mesh &sf_mesh,
                                            const AABBWrapper &aabb_wrapper,
                                            std::vector<REI> &rt_es_info) {
  // fetch all sf_mesh facets whose bbox has intersection with
  // the bbox of the dual polygon
  igl::parallel_for(
      rt_es_info.size(),
      [&](int i) {
        auto &rei = rt_es_info[i];
        if (rei.is_dual_poly_intersect) return;
        rei.poly_bbox_isect_triangles =
            aabb_wrapper.sf_polygon_box_intersections(sf_mesh,
                                                      rei.dual_polygon);
      },
      10000);

  // check if dual polygon has intersection with sf_mesh
  // update rei.is_dual_poly_intersect
  igl::parallel_for(rt_es_info.size(), [&rt_es_info, &sf_mesh](int i) {
    auto &rei = rt_es_info[i];
    if (rei.is_dual_poly_intersect || rei.poly_bbox_isect_triangles.size() == 0)
      return;
    for (const auto &f : rei.poly_bbox_isect_triangles) {
      // surface triangle
      std::vector<GEO::vec3> sf_triangle;
      for (int lv = 0; lv < sf_mesh.facets.nb_vertices(f); lv++)
        sf_triangle.push_back(
            sf_mesh.vertices.point(sf_mesh.facets.vertex(f, lv)));
      if (sf_triangle.size() != 3) log_and_throw("Sorry, triangle mesh only");

      // check intersection
      // triangulize the rei.dual_polygon
      GEO::vector<GEO::TriangleIsect> _;
      for (unsigned k = 1; k < rei.dual_polygon.size() - 1; k++) {
        // triangle (p1, p2, p3)
        GEO::vec3 p0 = to_vec3(rei.dual_polygon[0]);
        GEO::vec3 p1 = to_vec3(rei.dual_polygon[k]);
        GEO::vec3 p2 = to_vec3(rei.dual_polygon[k + 1]);
        bool is_intersect = false;
        try {
          is_intersect = GEO::triangles_intersections(
              p0, p1, p2, sf_triangle[0], sf_triangle[1], sf_triangle[2], _);
        } catch (const std::runtime_error &) {
          // some polygon triangle might be degenerated, do nothing
          break;
        }
        if (is_intersect) {
          rei.is_dual_poly_intersect = true;
          break;
        }
      }  // for poly_tris
      if (rei.is_dual_poly_intersect) break;
    }  // for rei.poly_bbox_isect_triangles
  });
}

// This function has been parallelized
void get_RT_edge_dual_info(const GEO::Mesh &sf_mesh,
                           const AABBWrapper &aabb_wrapper,
                           RegularTriangulationNN &rt) {
  logger().debug("[RT dual] fetching dual info of RT edges ...");
  int eid = 0;
  // edge is given by a cell and two indices (the edge (i,j) of a cell c is the
  // edge whose endpoints are the vertices of c with indices i and j).
  for (Finite_edges_iterator_rt fet = rt.finite_edges_begin();
       fet != rt.finite_edges_end(); fet++) {
    REI rei;
    // fetch edge 2 vertices info
    rei.valid_vs_tags[0] = fet->first->vertex(fet->second)->info().tag;
    rei.valid_vs_tags[1] = fet->first->vertex(fet->third)->info().tag;
    rei.all_vs_tags[0] = fet->first->vertex(fet->second)->info().all_tag;
    rei.all_vs_tags[1] = fet->first->vertex(fet->third)->info().all_tag;
    std::sort(rei.valid_vs_tags.begin(), rei.valid_vs_tags.end());
    std::sort(rei.all_vs_tags.begin(), rei.all_vs_tags.end());

    // fetch dual polygon info and check inifinity
    // incident_cells - circulator over
    // all cells incident to a given edge
    bool all_finite_inside = true;
    Cell_circulator_rt cc = rt.incident_cells(*fet);
    do {
      int tid = cc->info().id;
      if (rt.is_infinite(cc)) {
        all_finite_inside = false;
        break;
      }
      const auto &tet_info = rt.rt_ts_info.at(tid);
      rei.adj_tets.push_back(tid);
      rei.dual_polygon.push_back(tet_info.center);
      rei.is_dual_poly_vs_outside.push_back(tet_info.is_dual_point_outside);
      if (!tet_info.is_dual_point_outside) {
        // if any dual vertex (of tet) is inside, then
        // dual polygon (of edge) must intersects with object
        rei.is_dual_poly_intersect = true;
      }
      cc++;
    } while (cc != rt.incident_cells(*fet));
    // do not save edge
    if (!all_finite_inside) continue;

    // skip more calculations of polygon & mesh intersections
    // for making calculation faster
    // 1. store rei.adj_facets
    // 2. update rei.is_dual_poly_intersect
    int num_adj_tets = rei.adj_tets.size();
    for (int c = 0; c < num_adj_tets; c++) {
      auto &tid = rei.adj_tets[c];
      auto &tid_next = rei.adj_tets[(c + 1) % num_adj_tets];
      const auto &adj_fs = rt.rt_ts_info.at(tid).adj_facets;
      const auto &adj_fs_next = rt.rt_ts_info.at(tid_next).adj_facets;
      std::vector<int> isect_fs;
      set_intersection(adj_fs, adj_fs_next, isect_fs);
      if (isect_fs.size() != 1)
        log_and_throw("Two tets cannot adjacent more than 1 facet");
      int fid = isect_fs[0];
      rei.adj_facets.push_back(fid);
      const auto &f_info = rt.rt_fs_info.at(fid);
      if (f_info.is_dual_seg_intersect) {
        // if any dual segment (of facet) intersects, then
        // dual polygon (of edge) must intersects with object
        rei.is_dual_poly_intersect = true;
      }
    }
    // save anyway
    rei.id = eid++;
    rt.rt_es_info.push_back(rei);
  }

  // check if dual polygon has intersection with sf_mesh
  // update REI::is_dual_poly_intersect
  get_RT_edge_dual_poly_sf_intersections(sf_mesh, aabb_wrapper, rt.rt_es_info);

  logger().debug("[RT dual] fetched dual info of RT edges: {}",
                 rt.rt_es_info.size());
}

// Here we need to store the dual info of RT
// mostly for RT face, we need to store the dual segment
// here we need to know if each dual vertex is inside/outside
// but winding number is not consistent with AABB wrapper
// so we use a `evenlop_eps` as a narrow band around surface
// and consider any vertex in within this evenlop_eps as outside
void generate_RT_dual_info(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,  // for winding filter inside/outside
    const GEO::Mesh &sf_mesh, const AABBWrapper &aabb_wrapper,
    RegularTriangulationNN &rt, bool is_debug) {
  logger().debug("[RT dual] calling generate_RT_dual_info ...");
  igl::Timer timer;
  timer.start();
  get_RT_tet_dual_info(V, F, aabb_wrapper, rt);
  logger().info("[RT dual] dual tet info took {}s",
                timer.getElapsedTimeInSec());
  get_RT_face_dual_info(aabb_wrapper, rt);
  logger().info("[RT dual] dual face info took {}s",
                timer.getElapsedTimeInSec());
  get_RT_edge_dual_info(sf_mesh, aabb_wrapper, rt);
  logger().info("[RT dual] dual edge info took {}s",
                timer.getElapsedTimeInSec());
  logger().info("[RT dual] dual_intersection took {}s",
                timer.getElapsedTimeInSec());
  if (is_debug) rt.print_dual_info();
}

void generate_DT_CGAL(const std::vector<double> &pts,
                      const std::vector<double> &pts_normals,
                      const std::set<int> &seed_is_deleted,
                      const std::set<int> &feature_points,
                      DelaunayTriangulation &dt,
                      std::vector<MVertex> &all_medial_spheres, bool is_debug) {
  if (is_debug)
    logger().debug("----- generating Delaunay Triangulation using CGAL");
  const int dim = 3;
  const int nb_pts = pts.size() / dim;
  // sanity check
  if (pts_normals.size() != pts.size())
    log_and_throw("ERROR: pts_normals.size() != nb_pts");

  //////////
  // ATTENTION:
  // some inserted weighted points points can be hidden and do not result in
  // vertices in the triangulation, (coverted by other weighted points)
  // rt.number_of_vertices() <= wpts.size()
  //////////
  dt.clean();
  //////////
  // nb_pts = seed_points + 8 bbx (includes deleted points)
  dt.set_largest_v_tag(nb_pts);
  for (int i = 0; i < nb_pts; i++) {
    // we handle seed_is_deleted here becuase we want
    // dt vertices tags map to seed_points
    if (seed_is_deleted.find(i) != seed_is_deleted.end()) continue;

    Point_dt p(pts[i * dim], pts[i * dim + 1], pts[i * dim + 2]);
    // logger().debug("add mat vertex ({},{},{})",
    // 	pts[i*dim], pts[i*dim+1], pts[i*dim+2]);
    Vertex_handle_dt vh;
    vh = dt.insert(p);
    if (vh == nullptr) continue;
    vh->info().tag =
        i;  // this should be mapping to seed_points if not calculating LFS
    vh->info().pos = Vector3(p[0], p[1], p[2]);
    vh->info().normal = Vector3(pts_normals[i * dim], pts_normals[i * dim + 1],
                                pts_normals[i * dim + 2]);

    // not for LFS, since we ignore feature points
    if (feature_points.find(i) != feature_points.end())
      vh->info().is_on_feature = true;
  }
  assert(dt.is_valid());
  assert(dt.dimension() == 3);
  dt.print_info();

  // Step 1: compute circumcenters
  // Step 2: store radius^2 as weight for circumcenters
  std::vector<int> v_tags(4, -1);
  std::vector<GEO::vec3> v_pos_vec3(4);
  std::ostringstream stringStream;
  all_medial_spheres.clear();
  for (Finite_cells_iterator_dt fci = dt.finite_cells_begin();
       fci != dt.finite_cells_end(); fci++) {
    int all_idx = all_medial_spheres.size();
    // this is not identical during each run
    // since finite_cells_begin() starts with arbitrary cell
    // use fci->info().identifier instead
    fci->info().tag = all_idx;

    for (int i = 0; i < 4; i++) {
      v_tags[i] = fci->vertex(i)->info().tag;
      v_pos_vec3[i] = to_vec3(fci->vertex(i)->info().pos);
    }

    ////////////
    // create unique identifier
    // "x-x-x-x", each x is the tag of vertice (DVI::tag), sorting from low to
    // high
    std::sort(v_tags.begin(), v_tags.end());
    for (int i = 0; i < 4; i++) {
      stringStream << v_tags[i];
      if (i != 3) stringStream << "-";
    }
    fci->info().identifier = stringStream.str();
    stringStream.str("");  // clear

    // // calculate circumcenter and radius
    // Point_dt cent = CGAL::circumcenter(dt.tetrahedron(fci));
    // fci->info().center = Vector3(cent[0], cent[1], cent[2]);
    // fci->info().sq_radius = CGAL::squared_distance(cent,
    // fci->vertex(0)->point());

    // use geogram
    // default GEO::Geom::tetra_circum_center() seems not happy
    GEO::vec3 cent;
    bool is_good_tet = true;
    // cent = GEO::Geom::tetra_circum_center(v_pos_vec3[0], v_pos_vec3[1],
    // 								v_pos_vec3[2],
    // v_pos_vec3[3]);
    is_good_tet = tetra_circum_center(v_pos_vec3[0], v_pos_vec3[1],
                                      v_pos_vec3[2], v_pos_vec3[3], cent);
    if (!is_good_tet) {
      if (is_debug) {
        logger().error(
            "Error while computing circumsphere for seeds: {}, skip...",
            v_tags);
        for (int i = 0; i < 4; i++)
          logger().error("seed {}: ({},{},{})", v_tags[i], v_pos_vec3[i][0],
                         v_pos_vec3[i][1], v_pos_vec3[i][2]);
      }
      continue;
    }
    fci->info().center = Vector3(cent[0], cent[1], cent[2]);
    fci->info().sq_radius = GEO::Geom::distance2(cent, v_pos_vec3[0]);

    MVertex nf_mat_p(all_idx, fci->info().center[0], fci->info().center[1],
                     fci->info().center[2], fci->info().sq_radius);
    nf_mat_p.identifier = fci->info().identifier;
    // std::copy(v_tags.begin(), v_tags.end(),
    // std::back_inserter(nf_mat_p.dt_seed_ids)); push back 4 tet indices &
    // positions
    for (int i = 0; i < 4; i++) {
      nf_mat_p.dt_seed_ids[i] = v_tags[i];
      nf_mat_p.dt_seed_pos[i] = fci->vertex(i)->info().pos;
      nf_mat_p.dt_seed_normals[i] = fci->vertex(i)->info().normal;
    }
    // once push_back, nf_mat_p is just a copy
    all_medial_spheres.push_back(nf_mat_p);
  }
  dt.set_largest_cell_tag(all_medial_spheres.size());
}

//////////////////////////////////////////////////////////// DtIO

// for debug, write valid ma spheres to file
void DtIO::save_valid_medial_spheres(
    const std::vector<MVertex> &all_medial_spheres,
    const std::vector<int> &valid_medial_spheres) {
  // format: center[0] center[1] center[2] sq_radius

  if (all_medial_spheres.empty() && valid_medial_spheres.empty()) return;
  std::string output_path = "valid_medial_spheres.mvertex";

  logger().debug("Writing seeds to {} ...", output_path);
  std::ofstream f(output_path);
  for (int valid_idx = 0; valid_idx < valid_medial_spheres.size();
       valid_idx++) {
    const int all_idx = valid_medial_spheres.at(valid_idx);
    const MVertex &mat_p = all_medial_spheres.at(all_idx);

    // if (!mat_p.is_on_s_edge)
    // 	continue;

    // f << mat_p.tag << " ";
    for (int j = 0; j < 3; j++) {
      f << mat_p.pos[j] << " ";
    }
    f << mat_p.sq_radius << " ";

    f << std::endl;
  }
  f.close();
}

void DtIO::read_valid_medial_spheres(
    std::vector<MVertex> &all_medial_spheres,
    std::vector<int> &valid_medial_spheres,
    std::map<int, int> &all_to_valid_medial_spheres) {
  std::string input_path = "valid_medial_spheres.mvertex";
  logger().debug("Loading seeds from {} ...", input_path);
  all_medial_spheres.clear();
  valid_medial_spheres.clear();

  std::ifstream f(input_path.c_str());

  int tag = 0;
  double p0, p1, p2, sq_radius;

  while (f >> p0 >> p1 >> p2 >> sq_radius) {
    // if (sq_radius >= 0.0001) {
    // 	// sq_radius -= sq_radius*0.02;
    // 	// p1 -= p1*0.005;
    // 	double radius = std::sqrt(sq_radius);
    // 	radius = radius*0.9;
    // 	sq_radius = std::pow(radius, 2);
    // }

    MVertex mat_p(tag, p0, p1, p2, sq_radius);
    all_medial_spheres.push_back(mat_p);
    valid_medial_spheres.push_back(tag);
    all_to_valid_medial_spheres[tag] = tag;
    tag++;
  }

  logger().debug("#all_medial_spheres: {}", all_medial_spheres.size());
  logger().debug("#valid_medial_spheres: {}", valid_medial_spheres.size());
}

void DtIO::save_seeds_random(const std::string &output_path,
                             const std::vector<double> &seed_points,
                             const std::vector<double> &seed_normals) {
  if (seed_points.empty() && seed_normals.empty()) return;
  if (seed_points.size() != seed_normals.size())
    log_and_throw("ERROR: seed_points.size() != seed_normals.size()");

  const int dim = 3;
  const int seed_size = seed_points.size() / dim;
  logger().debug("Writing seeds to {} ...", output_path);
  std::ofstream f(output_path);
  // seed size
  // f -> non-feature seed points in 3D, normal in 3D
  // n -> feature seed points in 3D, normal in 3D
  f << seed_size << std::endl;
  for (int i = 0; i < seed_size; i++) {
    for (int j = 0; j < dim; j++) {
      f << seed_points[i * dim + j] << " ";
    }
    for (int j = 0; j < dim; j++) {
      f << seed_normals[i * dim + j] << " ";
    }
    f << std::endl;
  }
  f.close();
}

void DtIO::load_seeds_random(const std::string &input_path,
                             std::vector<double> &seed_points,
                             std::vector<double> &seed_normals) {
  logger().debug("Loading seeds from {} ...", input_path);
  seed_points.clear();
  seed_normals.clear();
  const int dim = 3;

  std::ifstream f(input_path.c_str());
  int nseeds;
  f >> nseeds;

  seed_points.resize(nseeds * dim);
  seed_normals.resize(nseeds * dim);

  for (int i = 0; i < nseeds; i++) {
    double x, y, z, n1, n2, n3;
    f >> x >> y >> z >> n1 >> n2 >> n3;
    seed_points[i * dim] = x;
    seed_points[i * dim + 1] = y;
    seed_points[i * dim + 2] = z;

    seed_normals[i * dim] = n1;
    seed_normals[i * dim + 1] = n2;
    seed_normals[i * dim + 2] = n3;
  }

  logger().debug("#seed_points: {}", seed_points.size());
  logger().debug("#seed_normals: {}", seed_normals.size());
}

void DtIO::save_seeds_after_CVT(const std::string &output_path,
                                const std::vector<double> &seed_points,
                                const std::vector<double> &seed_normals,
                                const std::set<int> &feature_points) {
  // save 3D point cloud
  bool is_save_seeds_only = true;
  if (is_save_seeds_only) {
    const int dim = 3;
    const int seed_size = seed_points.size() / dim;
    logger().debug("Writing seeds to {} ...", output_path);
    std::ofstream f(output_path);
    // x y z
    for (int i = 0; i < seed_size; i++) {
      for (int j = 0; j < dim; j++) {
        f << seed_points[i * dim + j] << " ";
      }
      f << std::endl;
    }
    f.close();
    return;
  }

  // deprecated
  if (seed_points.empty() && seed_normals.empty()) return;
  if (feature_points.empty())
    logger().debug("There is no feature points to save.");
  if (seed_points.size() != seed_normals.size())
    log_and_throw("ERROR: seed_points.size() != seed_normals.size()");

  const int dim = 3;
  const int seed_size = seed_points.size() / dim;
  logger().debug("Writing seeds to {} ...", output_path);
  std::ofstream f(output_path);
  // seed size
  // f -> non-feature seed points in 3D, normal in 3D
  // n -> feature seed points in 3D, normal in 3D
  f << seed_size << std::endl;
  for (int i = 0; i < seed_size; i++) {
    if (feature_points.find(i) != feature_points.end())
      f << "f ";
    else
      f << "n ";
    for (int j = 0; j < dim; j++) {
      f << seed_points[i * dim + j] << " ";
    }
    for (int j = 0; j < dim; j++) {
      f << seed_normals[i * dim + j] << " ";
    }
    f << std::endl;
  }
  f.close();
}

void DtIO::load_seeds(const std::string &input_path,
                      std::vector<double> &seed_points,
                      std::vector<double> &seed_normals,
                      std::set<int> &feature_points) {
  // deprecated!!!
  logger().debug("Loading seeds from {} ...", input_path);
  seed_points.clear();
  seed_normals.clear();
  feature_points.clear();
  const int dim = 3;

  std::ifstream f(input_path.c_str());
  int nseeds;
  f >> nseeds;

  seed_points.resize(nseeds * dim);
  seed_normals.resize(nseeds * dim);

  // feature
  for (int i = 0; i < nseeds; i++) {
    char ch;
    double x, y, z, n1, n2, n3;
    f >> ch >> x >> y >> z >> n1 >> n2 >> n3;

    if (ch == 'f') {
      feature_points.insert(i);
    }

    seed_points[i * dim] = x;
    seed_points[i * dim + 1] = y;
    seed_points[i * dim + 2] = z;

    seed_normals[i * dim] = n1;
    seed_normals[i * dim + 1] = n2;
    seed_normals[i * dim + 2] = n3;
  }

  logger().debug("#feature_points: {}", feature_points.size());
  logger().debug("#seed_points: {}", seed_points.size());
  logger().debug("#seed_normals: {}", seed_normals.size());
}

void DtIO::save_seeds_to_pts(const std::string &output_name,
                             const Parameters &params,
                             const std::vector<double> &seed_points,
                             bool is_normalize) {
  // save downsample as well
  double dc = params.downsample_percentage;
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << dc;
  std::string s_dc = stream.str();
  std::string output_path = "../out/pts/seeds_" + output_name + "_" + s_dc +
                            "_" + get_timestamp() + ".pts";

  logger().debug("Writing seeds to {} ...", output_path);
  ofstream f(output_path);
  if (!f.is_open()) {
    logger().critical("IOError: Could not open {}.", output_path);
    return;
  }

  const int seed_size = seed_points.size() / 3;

  if (is_normalize) {
    // normalization
    // min corner of AABB
    double xmin = numeric_limits<double>::max();
    double ymin = numeric_limits<double>::max();
    double zmin = numeric_limits<double>::max();
    // max corner of AABB
    double xmax = numeric_limits<double>::min();
    double ymax = numeric_limits<double>::min();
    double zmax = numeric_limits<double>::min();

    for (int i = 0; i < seed_size; i++) {
      const double &x = seed_points[i * 3];
      const double &y = seed_points[i * 3 + 1];
      const double &z = seed_points[i * 3 + 2];

      xmin = x < xmin ? x : xmin;
      ymin = y < ymin ? y : ymin;
      zmin = z < zmin ? z : zmin;

      xmax = x > xmax ? x : xmax;
      ymax = y > ymax ? y : ymax;
      zmax = z > zmax ? z : zmax;
    }

    logger().debug("Rescaling & Centering");
    logger().debug("min corner:[{},{},{}]", xmin, ymin, zmin);
    logger().debug("max corner:[{},{},{}]", xmax, ymax, zmax);

    double size = max(xmax - xmin, max(ymax - ymin, zmax - zmin)) / 10;

    double xcenter = (xmax + xmin) * 0.5;
    double ycenter = (ymax + ymin) * 0.5;
    double zcenter = (zmax + zmin) * 0.5;

    for (int i = 0; i < seed_size; i++) {
      f << (seed_points[i * 3] - xcenter) / size << " "
        << (seed_points[i * 3 + 1] - ycenter) / size << " "
        << (seed_points[i * 3 + 2] - zcenter) / size << endl;
    }
  } else {
    // no normalization
    for (int i = 0; i < seed_size; i++) {
      f << seed_points[i * 3] << " " << seed_points[i * 3 + 1] << " "
        << seed_points[i * 3 + 2] << endl;
    }
  }

  f.close();
  logger().debug("Done.");
  return;
}

// same as save_seeds_to_pts(), just in .ma format
void DtIO::save_seeds_to_ma(const std::string &output_name,
                            const Parameters &params,
                            const std::vector<double> &seed_points) {
  // save downsample as well
  double dc = params.downsample_percentage;
  std::stringstream stream;
  stream << std::fixed << std::setprecision(2) << dc;
  std::string s_dc = stream.str();
  std::string output_path = "../out/pts/seeds_" + output_name + "_" + s_dc +
                            "_" + get_timestamp() + ".ma";

  logger().debug("Writing seeds to {} ...", output_path);
  ofstream fout(output_path);
  if (!fout.is_open()) {
    logger().critical("IOError: Could not open {}.", output_path);
    return;
  }

  const int seed_size = seed_points.size() / 3;
  fout << seed_size << " " << 0 << " " << 0 << std::endl;

  for (int i = 0; i < seed_size; i++) {
    fout << "v " << seed_points[i * 3] << " " << seed_points[i * 3 + 1] << " "
         << seed_points[i * 3 + 2] << " " << 0 << endl;
  }
  fout.close();
  logger().debug("Done.");
  return;
}

void DtIO::save_all_mat_vertices(
    const std::string &output_path,
    const std::vector<double> &all_weighted_mat_pts_cgal,
    const std::vector<bool> &is_locking) {
  if (all_weighted_mat_pts_cgal.empty() && is_locking.empty()) return;

  int dim = 4;
  int all_size = all_weighted_mat_pts_cgal.size() / dim;
  if (all_size != is_locking.size())
    log_and_throw("ERROR: wrong during save_all_mat_vertices");

  logger().debug("Writing all MAT vertices to {} ...", output_path);
  std::ofstream f(output_path);
  // all_size
  // is_lock, x, y, z, weight
  f << all_size << std::endl;
  for (int i = 0; i < all_size; i++) {
    if (is_locking[i])
      f << 1 << " ";
    else
      f << 0 << " ";
    for (int j = 0; j < dim; j++) {
      f << all_weighted_mat_pts_cgal[i * dim + j] << " ";
    }
    f << std::endl;
  }
}

void DtIO::load_all_mat_vertices(const std::string &input_path,
                                 std::vector<double> &all_weighted_mat_pts_cgal,
                                 std::vector<bool> &is_locking) {
  logger().debug("Loading all MAT vertices from {} ...", input_path);
  all_weighted_mat_pts_cgal.clear();
  is_locking.clear();
  int dim = 4;

  std::ifstream f(input_path.c_str());
  int ntotal;
  f >> ntotal;

  all_weighted_mat_pts_cgal.resize(ntotal * dim);
  is_locking.resize(ntotal);

  for (int i = 0; i < ntotal; i++) {
    double ch, x, y, z, r;
    f >> ch >> x >> y >> z >> r;
    is_locking[i] = ch == 1. ? true : false;
    all_weighted_mat_pts_cgal[i * dim] = x;
    all_weighted_mat_pts_cgal[i * dim + 1] = y;
    all_weighted_mat_pts_cgal[i * dim + 2] = z;
    all_weighted_mat_pts_cgal[i * dim + 3] = r;
  }
  // logger().debug("ntotal: {}", ntotal);
  logger().debug("is_locking: {}", is_locking.size());
  logger().debug("all_weighted_mat_pts_cgal: {}",
                 all_weighted_mat_pts_cgal.size());
}

}  // namespace matfp