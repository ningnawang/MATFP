// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/LFS.h"

#include <igl/Timer.h>

#include "matfp/Triangulation.h"
#include "matfp/WindingFilter.h"
#include "matfp/geogram/mesh/mesh_sampling.h"

namespace matfp {

void get_bb_corners(Parameters& params, const std::vector<Vector3>& vertices,
                    Vector3& min, Vector3& max) {
  min = vertices[0];
  max = min;

  for (size_t j = 0; j < vertices.size(); j++) {
    for (int i = 0; i < 3; i++) {
      min(i) = std::min(min(i), vertices[j][i]);
      max(i) = std::max(max(i), vertices[j][i]);
    }
  }

  const Scalar dis = std::max(params.init_edge_length, params.eps_input * 2);
  for (int j = 0; j < 3; j++) {
    min[j] -= dis;
    max[j] += dis;
  }

  logger().debug("min = {} {} {}", min[0], min[1], min[2]);
  logger().debug("max = {} {} {}", max[0], max[1], max[2]);

  // params.bbox_min = min;
  // params.bbox_max = max;

  // params.bbox_min = Vector3(min[0], min[1], min[2]);
  // params.bbox_max = Vector3(max[0], max[1], max[2]);
}

// mainly for updating:
// - lfs_kd_tree
// - lfs_seeds
// - lfs_min_values
// - lfs_rho_values
void store_lfs_and_construct_kdtree(
    const std::vector<MVertex>& all_medial_spheres,
    const std::vector<double>& dt_seeds, const AABBWrapper& aabb_wrapper,
    GEO::NearestNeighborSearch_var& lfs_kd_tree, std::vector<double>& lfs_seeds,
    std::vector<double>& lfs_min_values, std::vector<double>& lfs_rho_values,
    const double r_sample, double& max_lfs_value, bool is_debug) {
  if (is_debug) logger().debug("calling store_lfs_and_construct_kdtree ...");
  lfs_seeds.clear();
  lfs_min_values.clear();
  std::map<int, int> dt_to_lfs_map;
  std::map<int, std::vector<double>>
      lfs_seed_idx_to_lfs_values;  // for selecting median
  for (const auto mat_p : all_medial_spheres) {
    // only store lfs on non-feature mat points
    if (mat_p.is_outside || mat_p.is_on_s_edge) continue;

    double minLFS = DBL_MAX;
    int lfs_seed_idx = -1;
    for (const auto dt_seed_idx : mat_p.dt_seed_ids) {
      // logger().debug("mat_p tag {} has dt_seed_idx: {}", mat_p.tag,
      // dt_seed_idx); logger().debug("------- mat_p tag {}", mat_p.tag);

      if (dt_to_lfs_map.find(dt_seed_idx) != dt_to_lfs_map.end()) {
        double radius = std::sqrt(mat_p.sq_radius);
        lfs_seed_idx = dt_to_lfs_map[dt_seed_idx];
        lfs_seed_idx_to_lfs_values[lfs_seed_idx].push_back(radius);
        // lfs_min_values[lfs_seed_idx] = std::min(
        // 	lfs_min_values[lfs_seed_idx],
        // 	radius
        // );
        // logger().debug("found dt_seed_idx: {} with lfs_seed_idx: {}, with
        // lfs_min_values[{}}]: {}", 	dt_seed_idx, lfs_seed_idx, lfs_seed_idx,
        // lfs_min_values[lfs_seed_idx]);
      } else {  // push new lfs seeds
        lfs_seed_idx = lfs_seeds.size() / 3;
        // logger().debug("not found dt_seed_idx: {}, add lfs_seed_idx: {}",
        // dt_seed_idx, lfs_seed_idx);
        dt_to_lfs_map[dt_seed_idx] = lfs_seed_idx;
        for (int i = 0; i < 3; i++)
          lfs_seeds.push_back(dt_seeds[dt_seed_idx * 3 + i]);

        double radius = std::sqrt(mat_p.sq_radius);

        // here we only store minimum radius size
        // minLFS = std::sqrt(mat_p.sq_radius);
        // we want to select median instead of minimum later
        lfs_seed_idx_to_lfs_values[lfs_seed_idx].push_back(radius);
        lfs_min_values.push_back(radius);
        // logger().debug("lfs_min_values[{}]: {}", lfs_min_values.size()-1,
        // lfs_min_values[lfs_min_values.size()-1]);
        // logger().debug("lfs_min_values[{}]: {} (lfs_seed_idx):",
        // lfs_seed_idx, lfs_min_values[lfs_seed_idx]);
        if (lfs_seed_idx != lfs_min_values.size() - 1) {
          logger().debug("lfs_min_values: {} != lfs_min_values.size(): {}",
                         lfs_seed_idx, lfs_min_values.size());
          log_and_throw("ERROR: lfs_seed_idx != lfs_min_values.size()-1");
        }
      }
      // logger().debug("dt_seed_idx {} has lfs_seed_idx {} has
      // lfs_min_values[lfs_seed_idx]: {}", 	dt_seed_idx, lfs_seed_idx,
      // lfs_min_values[lfs_seed_idx]);

      // // store maximum LFS for drawing colormap
      // max_lfs_value = std::max(max_lfs_value, lfs_min_values[lfs_seed_idx]);

    }  // for mat_p.dt_seed_ids done
  }    // for all_medial_spheres done

  // sanity check
  if (lfs_min_values.size() != lfs_seeds.size() / 3) {
    log_and_throw("ERROR: lfs_min_values.size() != lfs_seeds.size()/3");
  }

  // update lfs value to median, instead of minimum
  // to reduce the impact of spikes
  for (int i = 0; i < lfs_min_values.size(); i++) {
    std::vector<double>& values = lfs_seed_idx_to_lfs_values[i];
    lfs_min_values[i] = get_median(values);
    max_lfs_value = std::max(max_lfs_value, lfs_min_values[i]);
  }

  // store density/rho values
  lfs_rho_values.clear();
  lfs_rho_values.resize(lfs_min_values.size());
  for (int i = 0; i < lfs_min_values.size(); i++) {
    lfs_rho_values[i] = 1. / std::pow(r_sample * lfs_min_values[i], 4);
  }

  lfs_kd_tree = GEO::NearestNeighborSearch::create(3);
  lfs_kd_tree->set_points(lfs_seeds.size() / 3, lfs_seeds.data());

  if (is_debug) logger().debug("done store_lfs_and_construct_kdtree");
}

void lfs_generate_mesh_samples(const GEO::Mesh& sf_mesh,
                               const AABBWrapper& aabb_wrapper,
                               const size_t nb_points,
                               std::vector<Vector3>& sf_seeds,
                               std::vector<double>& points_new, bool is_debug) {
  if (is_debug) logger().debug("calling lfs_generate_mesh_samples...");
  sf_seeds.clear();
  const int dim = 3;
  int nb_new_points = nb_points;
  logger().info("Ideal #samples: {}", nb_points);
  if (nb_new_points < sf_mesh.vertices.nb()) {
    nb_new_points = sf_mesh.vertices.nb();
    if (is_debug)
      logger().info(
          "Current model meets our desired density, keep #nb_new_points {} ",
          nb_new_points);
  }

  //////////
  // Insert random samples
  std::vector<double> points_new_normals;  // no use
  GEO::Attribute<double> _;
  points_new.resize(nb_new_points * dim);
  points_new_normals.resize(nb_new_points * dim);
  // store new points and normals
  mesh_generate_random_samples_on_surface_with_normals<dim>(
      sf_mesh, points_new.data(), points_new_normals.data(), nb_new_points, _);
  if (is_debug) {
    logger().info("added {} random samples", nb_new_points);
    logger().debug("points_new.size(): {}", points_new.size() / 3);
    logger().debug("points_new_normals.size(): {}",
                   points_new_normals.size() / 3);
  }

  int removed_cnt = 0;
  for (int i = 0; i < points_new.size() / 3; i++) {
    // // delete seeds if too close to features
    // Vector3 p(points_new[i * 3], points_new[i * 3 + 1], points_new[i * 3 +
    // 2]); Scalar sq_dist = aabb_wrapper.get_sq_dist_to_s(p); if (sq_dist <=
    // 0.01) {
    //   removed_cnt++;
    //   continue;
    // }
    // store seeds and seeds normals
    sf_seeds.push_back(Vector3(points_new[i * 3], points_new[i * 3 + 1],
                               points_new[i * 3 + 2]));
  }
  if (is_debug) logger().debug("#all points size {}", sf_seeds.size());
}

void lfs_generate_dt(ThreeDimensionalShape& shape3D, bool is_debug) {
  igl::Timer timer;
  const int dim = 3;

  const std::vector<Vector3>& sf_seeds = shape3D.sf_seeds;
  const Eigen::MatrixXd& VI = shape3D.sf_mesh_wrapper.VI;
  const Eigen::MatrixXi& FI = shape3D.sf_mesh_wrapper.FI;

  std::vector<MVertex>& all_medial_spheres = shape3D.all_medial_spheres;
  DelaunayTriangulation& dt = shape3D.dt;
  std::vector<double> dt_vs, dt_vs_normals;

  bool& is_bbox_stored = shape3D.params.is_bbox_stored;
  std::vector<Vector3>& bb_points = shape3D.params.bb_points;

  // I: initial Delaunay vertices
  // step 1: add seeds
  int nb_pts = sf_seeds.size() + 8;
  for (int i = 0; i < sf_seeds.size(); i++) {
    for (int j = 0; j < dim; j++) {
      dt_vs.push_back(sf_seeds[i][j]);
      dt_vs_normals.push_back(0.);
    }
  }

  // step 1: add 8 bbox
  if (!is_bbox_stored) {
    Vector3 min, max;
    get_bb_corners(shape3D.params, sf_seeds, min, max);
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

  // step 3: calculate weight for non-feature mat points
  if (is_debug)
    logger().info("Initial DT points: {} = total {} + bbox 8", dt_vs.size() / 3,
                  sf_seeds.size());
  if (dt_vs.size() != dt_vs_normals.size())
    log_and_throw("ERROR: dt_vs.size() != dt_vs_normals.size()");

  // would store non-feature mat points
  std::set<int> _, _1;
  matfp::generate_DT_CGAL(dt_vs, dt_vs_normals, _, _1, dt, all_medial_spheres,
                          is_debug);

  timer.start();
  // step 4: in/out filter non-feature mat points
  std::vector<double> winding_num;
  std::vector<bool> is_outside;
  matfp::inout_filter_raw(all_medial_spheres, VI, FI, winding_num, is_outside,
                          dt);
  if (is_debug)
    logger().info("winding filter took {}s", timer.getElapsedTimeInSec());

  // step 6: add offset to inside spheres' radii for avoiding degeneracy
  // adding here instead of generate_RT_CGAL_MVertex()
  // is because function prepare_all_medial_spheres()
  // takes updated radius into account when deleting spheres
  for (auto& mat_p : all_medial_spheres) {
    if (mat_p.is_deleted || mat_p.is_outside) {
      continue;
    }
    mat_p.dilate_sphere_radius();
  }

  store_lfs_and_construct_kdtree(
      all_medial_spheres, dt_vs, shape3D.aabb_wrapper, shape3D.lfs_kd_tree,
      shape3D.lfs_seeds, shape3D.lfs_min_values, shape3D.rho_max_values,
      shape3D.params.r_sample, shape3D.max_lfs_value, is_debug);
}

////////////////////////////////////////////
// Main functions
////////////////////////////////////////////

void generate_lfs(ThreeDimensionalShape& shape3D, bool is_debug) {
  if (is_debug) logger().debug("calling generate_lfs...");
  lfs_generate_mesh_samples(shape3D.sf_mesh, shape3D.aabb_wrapper,
                            shape3D.params.init_nb_samples, shape3D.sf_seeds,
                            shape3D.init_seed_points, is_debug /*is_debug*/);
  lfs_generate_dt(shape3D, is_debug);
  if (is_debug) logger().debug("done generate_lfs");
}

}  // namespace matfp