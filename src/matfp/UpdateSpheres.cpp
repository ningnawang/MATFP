// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/UpdateSpheres.h"

#include "matfp/InternalFeatureAddition.h"
#include "matfp/IterateSpheres.h"

namespace matfp {

///////////
// check normal variance for all tangent planes
// if too small, then corresponding medial sphere
// is highly likely to be a spike
bool is_sphere_has_small_normal_variance(const MVertex& mat_p) {
  // compute centroid
  std::vector<Vector3> normals;
  Vector3 centroid(0., 0., 0.);
  for (auto& one_tan_pl : mat_p.tan_planes) {
    normals.push_back(one_tan_pl.normal);
    centroid += one_tan_pl.normal;
  }
  for (auto& one_cc_line : mat_p.tan_cc_lines) {  // also check tan_cc_lines
    normals.push_back(one_cc_line.normal);
    centroid += one_cc_line.normal;
  }
  centroid /= normals.size();
  centroid.normalize();  // project to unit sphere

  // here the spherical distance of two unit vectors d(p,q) = arcos(p.dot(q))
  double variance = 0.0;
  for (const auto& n : normals) {
    // the he domain of the arccosine function is from [−1, +1]
    // and the range is from [0, pi] radians (or from 0° to 180°).
    // double dist = std::acos(n.dot(centroid));
    double dist = spherical_dist_two_unit_vectors(n, centroid);
    double sq_dist = std::pow(dist, 2);
    variance += sq_dist;
  }
  variance /= normals.size();

  if (variance < 0.1)
    return true;
  else
    return false;
}

bool is_cc_part_has_small_normal_variance(
    const std::vector<Vector3> ref_fs_normals,
    const ConnectedComponent& one_cc_part) {
  const std::set<int>& ref_fs = one_cc_part.rpd_ref_fs_indices;
  // compute centroid
  std::vector<Vector3> normals;
  Vector3 centroid(0., 0., 0.);
  for (const auto& ref_f : ref_fs) {
    normals.push_back(ref_fs_normals[ref_f]);
    centroid += ref_fs_normals[ref_f];
  }
  centroid /= normals.size();
  centroid.normalize();  // project to unit sphere

  // here the spherical distance of two unit vectors d(p,q) = arcos(p.dot(q))
  double variance = 0.0;
  for (const auto& n : normals) {
    // the he domain of the arccosine function is from [−1, +1]
    // and the range is from [0, pi] radians (or from 0° to 180°).
    // double dist = std::acos(n.dot(centroid));
    double dist = spherical_dist_two_unit_vectors(n, centroid);
    double sq_dist = std::pow(dist, 2);
    variance += sq_dist;
  }
  variance /= normals.size();

  if (variance < 0.1)
    return true;
  else
    return false;
}

///////
// make sure all tan_elements has only one tangent point
// also delete those tangent planes that are marked as delete
// also remove sphere if contains no tangent planes
void purge_tan_planes(std::vector<MVertex>& all_medial_spheres) {
  logger().debug("start purging tangent planes for all spheres ...");
  std::set<int> purge_spikes;
  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) continue;

    if (mat_p.is_on_s_edge || mat_p.is_on_cc_edge) continue;

    if (mat_p.is_on_corner || mat_p.is_added_for_corner) continue;

    std::vector<TangentPlane>& tan_planes = mat_p.tan_planes;
    std::vector<TangentPlane> new_tan_planes;
    for (auto& tan_pl : tan_planes) {
      if (tan_pl.is_deleted)
        continue;
      else {
        tan_pl.update_points_with_centroid();
        tan_pl.sanity_check();
        new_tan_planes.push_back(tan_pl);
      }
    }
    tan_planes.clear();
    tan_planes = new_tan_planes;

    ////////////////////////////////////////////////
    // Purge Spikes
    //
    // delete if only contains 1 tangent plane
    // and not containing any concave line
    //
    // (WRONG?????)
    // Note: here we do not use size of cc_parts because
    //       spheres may not yet been pushed to ideal position
    //       so the collected cc_parts are not that reliable
    //       eg. some spheres are partually covered by another
    if (tan_planes.size() < 2 && mat_p.tan_cc_lines.size() < 1) {
      purge_spikes.insert(mat_p.tag);
      mat_p.is_deleted = true;
      continue;
    }

    // NOTE: this is too aggressive,
    //       many good spheres will be deleted during update
    //       because they might be covered by other bad sphere temperately
    // // this is a special case for the concave model
    // // TODO HERE
    // if (mat_p.cc_parts.size() < 1 && mat_p.tan_cc_lines.size() < 1) {
    // 	purge_spikes.insert(mat_p.tag);
    // 	mat_p.is_deleted = true;
    // 	continue;
    // }

    ///////////
    // check normal variance for all tangent planes
    // if too small, then corresponding medial sphere
    // is highly likely to be a spike
    if (is_sphere_has_small_normal_variance(mat_p)) {
      // logger().info("[Purge] mat_p {} has {} cc_parts, {} normals, centroid
      // ({},{},{}), variance: {}, DELETE", 	mat_p.tag,
      // mat_p.cc_parts.size(), 	normals.size(), centroid[0],
      // centroid[1], centroid[2], 	variance
      // );
      purge_spikes.insert(mat_p.tag);
      mat_p.is_deleted = true;
    }
  }  // for all_medial_spheres

  // if (purge_spikes.find(16162) != purge_spikes.end()) {
  // 	logger().error("mat {} has been purged!!!", 16162);
  // }

  logger().debug("[Purge] removed {} medial spheres as spikes",
                 purge_spikes.size());
}

// remove low quality spheres
// 1. if contains only < 2 cc_part and no concave line, and normal variance is
// small
// 2. for non-feature sphere, all its cc_parts are around concave lines
// 3. extrude too much from any tangent plane
void clean_spheres(const std::vector<Vector3> ref_fs_normals,
                   std::vector<MVertex>& all_medial_spheres,
                   bool is_clean_extrude, bool is_clean_no_cc,
                   bool is_remove_dup) {
  std::set<int> cleaned_spheres;
  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside || mat_p.is_on_bbox) continue;

    if (mat_p.is_on_s_edge || mat_p.is_on_cc_edge) continue;

    if (mat_p.is_on_corner || mat_p.is_added_for_corner) continue;

    if (mat_p.tan_planes.size() < 2 && mat_p.tan_cc_lines.size() < 1) {
      cleaned_spheres.insert(mat_p.tag);
      mat_p.is_deleted = true;
      continue;
    }

    // // some non-feature spheres are too close to features
    // if (mat_p.tan_planes.size() == 2 && mat_p.cc_parts.size() > 2) {
    //   cleaned_spheres.insert(mat_p.tag);
    //   mat_p.is_deleted = true;
    //   mat_p.is_deleted_int = DeletionType::D_CLEAN;
    //   continue;
    // }

    // clean sphere if it has no cc so type is unknown
    if (is_clean_no_cc) {
      if (mat_p.cc_parts.size() == 0 || mat_p.type == SphereType::T_UNK) {
        mat_p.is_deleted = true;
      }
    }

    // // type 1:
    // // if contains only < 2 cc_part and no concave line, and normal variance
    // is small
    // // then we clean it up
    // if (mat_p.cc_parts.empty()) {
    // 	logger().debug("[Clean] mat {} contains 0 cc_parts, delete", mat_p.tag);
    // 	cleaned_spheres.insert(mat_p.tag);
    // 	mat_p.is_deleted = true;
    // 	continue;
    // }
    // if (mat_p.cc_parts.size() == 1
    // 	&& mat_p.tan_cc_lines.empty()
    // 	&& is_cc_part_has_small_normal_variance(ref_fs_normals,
    // mat_p.cc_parts[0]) ) { 	logger().debug("[Clean] mat {} contains 1
    // cc_parts with small normal variance, delete", mat_p.tag);
    // 	cleaned_spheres.insert(mat_p.tag);
    // 	mat_p.is_deleted = true;
    // 	continue;
    // }

    // // type 2:
    // // if all cc_parts are around concave lines, for non-feature spheres
    // if (!mat_p.is_a_feature_sphere()) {
    // 	bool is_delete = true;
    // 	for (const auto& one_cc_part: mat_p.cc_parts) {
    // 		if (!one_cc_part.is_on_cc_line) {
    // 			is_delete = false;
    // 		}
    // 	}
    // 	if (is_delete) {
    // 		logger().debug("[Clean] mat {} only has RPD around concave lines
    // {},
    // delete", mat_p.tag); 		cleaned_spheres.insert(mat_p.tag);
    // mat_p.is_deleted = true; 		continue;
    // 	}
    // }

    // type 3:
    // check all tangent planes, if extrude too much, then delete
    // only check non-feature spheres
    if (is_clean_extrude) {
      // if (mat_p.is_on_intf_edge) continue;
      double extrude_threshold = 0.05;
      const Vector3& center = mat_p.pos;
      const double radius = std::sqrt(mat_p.sq_radius);
      for (const int& tan_pl_id : mat_p.selected_tan_pls) {
        auto& tan_pl = mat_p.tan_planes[tan_pl_id];
        tan_pl.update_extrude_ratio(mat_p.pos, mat_p.sq_radius);
        if (tan_pl.extrude_ratio > extrude_threshold) {
          logger().debug("[Clean] mat {} extrude tan_pl too much {}, delete",
                         mat_p.tag, tan_pl.extrude_ratio);
          cleaned_spheres.insert(mat_p.tag);
          mat_p.is_deleted = true;
          break;
        }
      }
    }

    // // check adjacent facets for tan_cc_lines as well
    // bool is_skip = false;
    // for (auto& tan_cc_line: mat_p.tan_cc_lines) {
    // 	const Vector3& point = tan_cc_line.tan_point;
    // 	for (const auto& adj_ref_f: tan_cc_line.all_adj_ref_fs) {
    // 		const Vector3& fnormal = ref_fs_normals[adj_ref_f];
    // 		double ratio = std::abs(radius - (point - center).dot(fnormal))
    // / radius; 		tan_cc_line.max_extrude_ratio =
    // std::max(tan_cc_line.max_extrude_ratio, ratio); 		if (ratio > 0.3)
    // { 			logger().debug("[Clean] mat {} extrude cc_line
    // too much
    // {}, delete", mat_p.tag, ratio); cleaned_spheres.insert(mat_p.tag);
    // mat_p.is_deleted = true; is_skip = true; break;
    // 		}
    // 		if (is_skip)
    // 			break;
    // 	}
    // 	if (is_skip)
    // 		break;
    // }

  }  // for all_medial_spheres
  logger().debug("[Clean] cleaned {} medial spheres: {}",
                 cleaned_spheres.size(), cleaned_spheres);

  // Remove duplicates
  if (is_remove_dup) {
    // if is_override = true: check all medial spehres, not only newly added
    std::vector<MVertex> sorted_partial_medial_spheres;
    int dup_cnt = remove_duplicated_medial_spheres(
        sorted_partial_medial_spheres, all_medial_spheres,
        true /*is_override*/);
  }
}

// no use
void purge_cc_parts(std::vector<MVertex>& all_medial_spheres) {
  logger().debug("[Purge] start purging medial spheres ...");
  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_outside) continue;

    // remove deleted cc parts
    std::vector<ConnectedComponent>& cc_parts = mat_p.cc_parts;
    std::vector<ConnectedComponent> new_cc_parts;
    for (auto& one_cc_part : cc_parts) {
      if (one_cc_part.is_deleted)
        continue;
      else {
        new_cc_parts.push_back(one_cc_part);
      }
    }  // for cc_parts
    cc_parts.clear();
    cc_parts = new_cc_parts;
  }  // for all_medial_spheres
}

void load_nearby_rpd_info_same_cc_part(
    const GEO::Mesh& rpd_mesh,
    const GEO::Attribute<GEO::index_t>& facet_region_attr,
    const std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_vs_bisectors,
    const std::vector<int>& valid_medial_spheres,
    const std::map<int, std::unordered_set<int>>& vertex_2_rpd_facets,
    const int& cur_f, const int& cur_seed_idx, ConnectedComponent& one_cc_part,
    std::vector<int>& cc_part_stack) {
  // loop to find other facets in the same cc part
  for (int c = rpd_mesh.facets.corners_begin(cur_f);
       c < rpd_mesh.facets.corners_end(cur_f); ++c) {
    int v = rpd_mesh.facet_corners.vertex(c);
    int v_next = rpd_mesh.facet_corners.vertex(
        rpd_mesh.facets.next_corner_around_facet(cur_f, c));
    // push rpd boundary vertices to given cc_part
    const std::set<GEO::index_t>& v_bisectors = rpd_vs_bisectors.at(v);
    if (v_bisectors.size() > 2) {
      one_cc_part.rpd_boundary_vertices.insert(v);
    }
    // found the rpd facet in the same connected components
    std::vector<int> common_f_ids;
    set_intersection(vertex_2_rpd_facets.at(v), vertex_2_rpd_facets.at(v_next),
                     common_f_ids);
    if (common_f_ids.size() != 2) {
      // logger().error("[RPD] We do not know how to handle non-manifold edge
      // for mat {}, facet {}, common_f_ids: {}", 	mat_p.tag, cur_f,
      // common_f_ids);
      continue;
    }
    for (const auto& adj_f : common_f_ids) {
      if (adj_f == cur_f) continue;
      const int adj_seed_idx = facet_region_attr[adj_f];
      if (adj_seed_idx == cur_seed_idx) {
        // cur_f and adj_f are in the same cc_part
        cc_part_stack.push_back(adj_f);
      } else {
        const int adj_tag = valid_medial_spheres.at(adj_seed_idx);
        one_cc_part.adj_sphere_tags.insert(adj_tag);
      }
    }
  }  // for lv in rpd_mesh.facets
}

// TODO: deprecate
void update_tan_cc_lines_using_RPD(
    const GEO::Mesh& rpd_mesh,
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres) {
  // if no cc line, just skip
  if (tan_cc_lines.empty()) return;

  const GEO::Attribute<GEO::index_t> facet_region_attr(
      rpd_mesh.facets.attributes(), "region");
  const GEO::Attribute<GEO::index_t> facet_ref_facet_attr(
      rpd_mesh.facets.attributes(),
      "ref_facet");  // for input mesh facet indices

  std::map<int, std::set<int>> valid_2_ref_facets;
  // for detecting center of concave line segment
  std::map<int, std::set<int>> vertex_2_ref_facets;

  // update valid_2_ref_facets [reuse]
  // here we store all reference facets, reuse later
  for (int f = 0; f < rpd_mesh.facets.nb(); f++) {
    const int seed_idx = facet_region_attr[f];  // = valid_idx
    const int ref_facet = facet_ref_facet_attr[f];
    valid_2_ref_facets[seed_idx].insert(ref_facet);
    for (int c = rpd_mesh.facets.corners_begin(f);
         c < rpd_mesh.facets.corners_end(f); ++c) {
      int v = rpd_mesh.facet_corners.vertex(c);
      vertex_2_ref_facets[v].insert(ref_facet);
    }
  }

  std::set<int> cc_related_spheres;  // for sanity check
  for (auto& pair : valid_2_ref_facets) {
    const int valid_idx = pair.first;
    std::set<int>& ref_facets = pair.second;
    MVertex& mat_p = all_medial_spheres.at(valid_medial_spheres.at(valid_idx));
    std::set<int> new_ref_faces;  // will replace ref_facets
    // mat_p.tan_cc_lines.clear();

    // tan_cc_lines size should be relative small (wrong!!!)
    // 2021-09-04 ninwang: tan_cc_lines size could be large
    for (const auto& one_cc_line : tan_cc_lines) {
      // check if cc_line already exist
      bool is_save_new_cc_line = true;
      for (auto& tan_cc_line : mat_p.tan_cc_lines) {
        if (tan_cc_line == one_cc_line) {
          is_save_new_cc_line = false;
          // one_cc_line.store_all_adj_ref_fs(new_ref_faces);
          // tan_cc_line.is_tan_point_updated = false;
          cc_related_spheres.insert(valid_idx);
          break;
        }
      }
      if (!is_save_new_cc_line) continue;

      // // to find
      // for (const auto& ref_fs_pair : one_cc_line.adj_ref_fs_pair) {
      //   const int f1 = ref_fs_pair[0];
      //   const int f2 = ref_fs_pair[1];
      //   if (ref_facets.find(f1) == ref_facets.end()) continue;
      //   if (ref_facets.find(f2) == ref_facets.end()) continue;
      //   // found new, store it
      //   mat_p.tan_cc_lines.push_back(one_cc_line);
      //   // mat_p.tan_cc_lines.back().is_tan_point_updated = false;
      //   mat_p.is_to_be_updated = true;
      //   cc_related_spheres.insert(valid_idx);
      //   // one_cc_line.store_all_adj_ref_fs(new_ref_faces);
      //   break;
      // }
    }
    ref_facets.clear();
    ref_facets = new_ref_faces;
  }  // for valid_2_ref_facets end

  // update tangent point X as center of concave segment
  std::set<int> cc_normal_updated_spheres;  // for sanity check
  for (int f = 0; f < rpd_mesh.facets.nb(); f++) {
    const int valid_idx = facet_region_attr[f];  // = valid_idx
    const int ref_facet = facet_ref_facet_attr[f];
    MVertex& mat_p = all_medial_spheres.at(valid_medial_spheres.at(valid_idx));
    if (mat_p.tan_cc_lines.empty()) continue;

    for (int c = rpd_mesh.facets.corners_begin(f);
         c < rpd_mesh.facets.corners_end(f); ++c) {
      int v1 = rpd_mesh.facet_corners.vertex(c);
      int v2 = rpd_mesh.facet_corners.vertex(
          rpd_mesh.facets.next_corner_around_facet(f, c));
      std::set<int> v1_ref_faces = vertex_2_ref_facets.at(v1);
      std::set<int> v2_ref_faces = vertex_2_ref_facets.at(v2);
      std::vector<int> inter_ref_faces;
      set_intersection<int>(v1_ref_faces, v2_ref_faces, inter_ref_faces);
      if (inter_ref_faces.size() != 2) continue;

      std::array<int, 2> v_ref_faces = {
          {inter_ref_faces[0], inter_ref_faces[1]}};
      std::sort(v_ref_faces.begin(), v_ref_faces.end());

      // // try to find
      // for (auto& one_cc_line : mat_p.tan_cc_lines) {
      //   for (const auto& ref_fs_pair : one_cc_line.adj_ref_fs_pair) {
      //     if (v_ref_faces == ref_fs_pair) {
      //       // update tangent concave segment
      //       const GEO::vec3& v1_p = rpd_mesh.vertices.point(v1);
      //       const GEO::vec3& v2_p = rpd_mesh.vertices.point(v2);
      //       one_cc_line.update_tangent_point(v1_p, v2_p);
      //       // one_cc_line.update_tangent_normal(); // update when creating a
      //       // new sphere
      //       one_cc_line.is_tan_point_updated = true;
      //       mat_p.is_to_be_updated = true;
      //       cc_normal_updated_spheres.insert(valid_idx);

      //       // sanity check
      //       // this sphere must exists in cc_related_spheres
      //       if (cc_related_spheres.find(valid_idx) ==
      //           cc_related_spheres.end()) {
      //         logger().error("cc_related_spheres: {}, not contain valid_idx
      //         {}",
      //                        cc_related_spheres, valid_idx);
      //         log_and_throw("ERROR");
      //       }
      //       break;
      //     }
      //   }
      // }  // for mat_p.tan_cc_lines
    }  // for rpd_mesh.facets.corners
  }    // for rpd_mesh.facets

  // // sanity check
  // // all spheres contains cc should be updated their cc normals
  // if (cc_normal_updated_spheres != cc_related_spheres) {
  // 	logger().error("cc_related_spheres: {}, cc_normal_updated_spheres: {}",
  // 		cc_related_spheres, cc_normal_updated_spheres);
  // 	log_and_throw("ERROR: cc_related_spheres != cc_normal_updated_spheres");
  // }
}

///////////////////////////////////////////////////////
// Main functions
// NOTE:
// 1. tangent elements are accumulated
// 2. cc_parts are not accumulated, will be cleared
//    and re-stored for each iteration with new RPD
//
// Note: Tangent planes of 1. feature sphere and 2. sphere added for corners
//       will be updated no matter the flag is_update_tan_pls
void update_tan_elements_and_cc_parts_using_RDP(
    const GEO::Mesh& rpd_mesh,
    const std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_vs_bisectors,
    const std::vector<Vector3> ref_fs_normals, const bool is_volumetric,
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_update_tan_pls) {
  const GEO::Attribute<GEO::index_t> facet_region_attr(
      rpd_mesh.facets.attributes(), "region");
  const GEO::Attribute<GEO::index_t> facet_ref_facet_attr(
      rpd_mesh.facets.attributes(),
      "ref_facet");  // for input mesh facet indices
  const GEO::Attribute<GEO::signed_index_t> facet_adj_region_attr(
      rpd_mesh.facets.attributes(), "adj_region");

  // Loop RPD mesh once to fetch all infos
  // for tan_cc_lines
  std::map<int, std::set<int>> valid_2_ref_facets;
  std::map<int, std::set<int>>
      vertex_2_ref_facets;  // for detecting center of concave line segment
  // for cc_parts
  std::map<int, std::unordered_set<int>> vertex_2_rpd_facets;
  for (int f = 0; f < rpd_mesh.facets.nb(); f++) {
    // for volumetric, we only care about surface
    if (is_volumetric) {
      const GEO::signed_index_t adj_idx = facet_adj_region_attr[f];
      if (adj_idx != -1) continue;
    }
    const int seed_idx = facet_region_attr[f];  // = valid_idx
    const int all_idx = valid_medial_spheres.at(seed_idx);
    const int ref_facet = facet_ref_facet_attr[f];
    MVertex& mat_p = all_medial_spheres.at(all_idx);

    valid_2_ref_facets[seed_idx].insert(ref_facet);
    for (int c = rpd_mesh.facets.corners_begin(f);
         c < rpd_mesh.facets.corners_end(f); ++c) {
      int v = rpd_mesh.facet_corners.vertex(c);
      vertex_2_ref_facets[v].insert(ref_facet);
      vertex_2_rpd_facets[v].insert(f);
    }
  }
  // Here we loop again because some sphere may not have
  // RPD, may covered by other sphere
  for (auto& mat_p : all_medial_spheres) {
    // NOTE: we refresh cc_part everytime for a new RPD
    //       cc_parts are not accumulated!!!
    mat_p.cc_parts.clear();
  }

  ////////////////////////////////
  // Step 2:
  // update tangent concave lines
  // [reuse] valid_2_ref_facets
  // valid_2_ref_facets will be updated,
  // only tan_cc_lines related reference facets will be stored.
  // used later to ignore those tangent planes
  // update_tan_cc_lines_using_RPD(rpd_mesh, tan_cc_lines, valid_medial_spheres,
  //                               all_medial_spheres);

  /////////////////////////////////////////////
  // Step 3:
  // 1. create new connected components
  //    cc_parts only store rpd adjacent info
  // 	  cc_parts will be cleared and updated every time
  // 2. crerate new tangent planes
  // 	  ignore those planes added to concave lines
  std::set<int> is_f_visited;
  std::vector<int> cc_part_stack;
  for (int f = 0; f < rpd_mesh.facets.nb(); f++) {
    // not reliable
    if (rpd_mesh.facets.nb_vertices(f) < 3) continue;
    if (is_f_visited.find(f) != is_f_visited.end()) continue;

    const int seed_idx = facet_region_attr[f];  // = valid_idx
    const int all_idx = valid_medial_spheres.at(seed_idx);
    MVertex& mat_p = all_medial_spheres.at(all_idx);
    int ref_facet;    // = facet_ref_facet_attr[f];
    Vector3 fnormal;  // = ref_fs_normals[ref_facet]; // use reference/input
                      // facets' normals
    ConnectedComponent one_cc_part;
    cc_part_stack.push_back(f);

    // if (all_idx == 1003) {
    // 	logger().debug("[CC_part] loading new cc_part to mat {}", all_idx);
    // }

    // In the while loop,
    // we will: 3.1. update tangent planes
    // 			3.2. update cc_part_stack
    while (!cc_part_stack.empty()) {
      int cur_f = cc_part_stack.back();
      cc_part_stack.pop_back();
      if (is_f_visited.find(cur_f) != is_f_visited.end()) continue;
      is_f_visited.insert(cur_f);

      if (rpd_mesh.facets.nb_vertices(cur_f) < 3) continue;

      ref_facet = facet_ref_facet_attr[cur_f];
      // computing RPD in parallel gives us -1 ref_facet
      // so we cannot update tangent planes
      if (ref_facet < 0 || ref_facet >= ref_fs_normals.size()) {
        // logger().error("ERROR: mat_p {} has RPD f {} with ref_facet {}, set
        // is_update_tan_pls = false", mat_p.tag, cur_f, ref_facet);
        is_update_tan_pls = false;
        // log_and_throw("ERROR");
      } else {
        fnormal =
            ref_fs_normals[ref_facet];  // use reference/input facets' normals
      }
      GEO::vec3 c = GEO::Geom::mesh_facet_center(rpd_mesh, cur_f);
      Vector3 fcentroid(c[0], c[1], c[2]);

      one_cc_part.rpd_facet_indices.insert(cur_f);
      one_cc_part.rpd_ref_fs_indices.insert(ref_facet);

      // // Check Concave Lines
      // bool is_erase_tan_pl_for_cc = false;
      // if (mat_p.tan_cc_lines.size() > 0) {
      //   for (const auto& tan_cc_line : mat_p.tan_cc_lines) {
      //     // Tangent Planes
      //     // Only checking two reference faces is not enough
      //     // some RPD facet is not on either reference faces
      //     // here we check the normal for tan_cc_lines
      //     if (tan_cc_line.is_normal_covered_by_adj_fs(fnormal)) {
      //       is_erase_tan_pl_for_cc = true;
      //     }
      //     // Connected Components
      //     // mark that cc_part that cross any concave line
      //     if (tan_cc_line.all_adj_ref_fs.find(ref_facet) !=
      //         tan_cc_line.all_adj_ref_fs.end()) {
      //       one_cc_part.is_on_cc_line = true;
      //     }
      //   }
      // }  // if tan_cc_lines.size() > 0

      ////////////////////////////////////////
      // 3.1 Tangent Plane: check and update
      //
      // any tangent planes that has same normal
      // as this plane for cc will be deleted
      // do not save if already saved for tan_cc_lines
      // loop existing tanget elements
      //
      if (mat_p.is_a_feature_sphere()) {
        is_update_tan_pls = false;
      }
      // Note: sphere added for corners whose
      //       tangent planes will be updated no matter flag is_update_tan_pls
      if (is_update_tan_pls || mat_p.is_added_for_corner) {
        bool is_save_new_tan_pl = true;
        if (mat_p.is_normal_covered_by_tan_pls(
                fnormal)) {  // do not save current
          is_save_new_tan_pl = false;
        }
        // for (auto& tan_pl : mat_p.tan_planes) {
        //   if (tan_pl.is_same_normal(fnormal)) {  // do not save current
        //     is_save_new_tan_pl = false;
        //     // if (is_erase_tan_pl_for_cc) {
        //     // 	tan_pl.is_deleted = true;
        //     // 	continue;
        //     // }
        //     // if same plane, update tangent point
        //     tan_pl.push_new_point(fcentroid);
        //     // mat_p.is_to_be_updated = true;
        //     break;
        //   }
        // }
        // create new tangent plane if different
        if (is_save_new_tan_pl) {
          mat_p.is_to_be_updated = true;
          TangentPlane tan_pl1(fnormal);
          tan_pl1.push_new_point(fcentroid);
          mat_p.tan_planes.push_back(tan_pl1);

          // if (mat_p.tag == 32971) {
          //     logger().error("------- mat {} adding tangent plane",
          //     mat_p.tag); tan_pl1.print_info();
          // }
        }
      }  // is_update_tan_pls

      //////////////////////////////////////
      // 3.2 Connected Component: update cc_part_stack
      load_nearby_rpd_info_same_cc_part(
          rpd_mesh, facet_region_attr, rpd_vs_bisectors, valid_medial_spheres,
          vertex_2_rpd_facets, cur_f, seed_idx, one_cc_part, cc_part_stack);

    }  // while cc_part_stack
    mat_p.cc_parts.push_back(one_cc_part);
  }  // for rpd_mesh facets

  ///////
  // 1. make sure all tan_elements has only one tangent point
  // 2. delete those tangent planes that are marked as delete
  // 3. delete sphere if normal variance is small
  // there is no need to purge cc_parts because we update it
  // once we recompute the RDP
  purge_tan_planes(all_medial_spheres);
}

bool is_sphere_need_to_be_updated(MVertex& mat_p) {
  double extrude_threshold = 0.05;
  bool is_update = false;
  // check all tangent planes, if extrude too much, then delete
  const Vector3& center = mat_p.pos;
  const double radius = std::sqrt(mat_p.sq_radius);
  for (auto& tan_pl : mat_p.tan_planes) {
    tan_pl.update_extrude_ratio(mat_p.pos, mat_p.sq_radius);
    if (!is_update && tan_pl.extrude_ratio > extrude_threshold) {
      is_update = true;
      logger().debug(
          "[IS_UPDATE] mat_p {} has tan_pl.extrude_ratio {}, going to be "
          "updated",
          mat_p.tag, tan_pl.extrude_ratio);
      // do not break here
      // we want to update all extrude ratio
    }
  }
  return is_update;
}

// will remove spheres that close to concave lines
void update_spheres(const GEO::Mesh& sf_mesh,
                    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
                    const AABBWrapper& aabb_wrapper,
                    std::vector<MVertex>& all_medial_spheres, bool is_debug) {
  // make sure all new updated spheres are unique
  std::vector<MVertex> sorted_partial_medial_spheres;
  std::set<int> updated_spheres, deleted_spheres;
  for (MVertex& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) continue;
    if (mat_p.is_on_s_edge || mat_p.is_on_corner) continue;
    // if (mat_p.is_on_cc_edge || mat_p.is_added_for_corner) continue;

    // only update T_3 or higher
    if (mat_p.tan_planes.size() <= 2) continue;

    ////////////
    // Do not delete when mat_p.tan_planes.empty(),
    // spheres added for corners
    // when mat_p.is_added_for_corner == true
    // do not have tangent planes
    if (!mat_p.is_to_be_updated && !is_sphere_need_to_be_updated(mat_p))
      continue;

    // do not update for spheres that are added only for concave lines
    // for those spheres, we use shrink_spheres() instead
    if (mat_p.type == SphereType::T_2_c || mat_p.type == SphereType::T_N_c)
      continue;

    mat_p.is_to_be_updated = false;
    mat_p.is_sphere_updated = true;
    updated_spheres.insert(mat_p.tag);
    if (is_debug)
      logger().debug("[Update] mat_p {} is going to be updated!!!", mat_p.tag);
    bool is_success = iterate_sphere(sf_mesh, ref_fs_pairs_not_cross,
                                     aabb_wrapper, mat_p, is_debug);

    if (!is_success) {
      if (is_debug)
        logger().error(
            "[Update] mat {} update failed, tan_pls {}, tan_cc_lines {}",
            mat_p.tag, mat_p.tan_planes.size(), mat_p.tan_cc_lines.size());
      mat_p.is_deleted = true;
      deleted_spheres.insert(mat_p.tag);
      continue;
    }

    // Post checking spheres
    // if the q point is too close to concave lines, then delete
    // this is because we have inserted enough spheres around concave lines
    // by setting pin point p on concave lines and shrink
    for (const auto& tan_pl : mat_p.tan_planes) {
      double p_sq_dist_2cc = aabb_wrapper.get_sq_dist_to_cc(tan_pl.points[0]);
      if (p_sq_dist_2cc <= SQ_DIST_TO_CC) {
        // logger().error("[Update] mat {} too close to concave lines, delete",
        //                mat_p.tag);
        mat_p.is_deleted = true;
        deleted_spheres.insert(mat_p.tag);
        break;
      }
    }  // for mat_p.tan_planes

    // All good
    // remember to update dilated radius
    mat_p.dilate_sphere_radius();
    if (mat_p.tan_planes.size() == 2) {
      mat_p.type = SphereType::T_2;
    } else {
      mat_p.type = SphereType::T_N;
    }

    if (is_debug) {
      logger().debug(
          "[Update] mat {}, type {}, new pos ({},{},{}) -> ({},{},{}), "
          "sq_radius "
          "{} -> "
          "{}, "
          "sq_radius_dilated: {}",
          mat_p.tag, mat_p.type, mat_p.old_pos[0], mat_p.old_pos[1],
          mat_p.old_pos[2], mat_p.pos[0], mat_p.pos[1], mat_p.pos[2],
          mat_p.old_sq_radius, mat_p.sq_radius, mat_p.sq_radius_dilated);
    }
    sorted_partial_medial_spheres.push_back(mat_p);

  }  // for all_medial_spheres

  logger().debug("[Update] updated {} spheres, deleted {} spheres",
                 updated_spheres.size(), deleted_spheres.size());
  // Remove duplicates
  // sort spheres and delete duplicated spheres
  remove_duplicated_medial_spheres(sorted_partial_medial_spheres,
                                   all_medial_spheres);

  purge_tan_planes(all_medial_spheres);
}

}  // namespace matfp