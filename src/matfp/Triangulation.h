// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/delaunay/delaunay.h>
#include <geogram/points/kd_tree.h>
#include <geogram/points/nn_search.h>

#include <limits>

#include "matfp/AABBWrapper.h"
#include "matfp/Common.h"
#include "matfp/Logger.h"
#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/Parameters.h"
#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/Types/RegularTriangulation.h"

namespace matfp {
// generate Delaunay Triangulation using CGAL
// output weighted points for Regular Triangulation (Weighted Delaunay)
void generate_DT_CGAL(const std::vector<double>& pts,
                      const std::vector<double>& pts_normals,
                      const std::set<int>& seed_is_deleted,
                      const std::set<int>& feature_points,
                      DelaunayTriangulation& dt,
                      std::vector<MVertex>& all_medial_spheres, bool is_debug);

// Regular triangulation
// generate RT for computing MAT using dual PD
void generate_RT_for_dual_PD(std::vector<MVertex>& all_medial_spheres,
                             std::vector<int>& valid_medial_spheres,
                             const std::vector<double>& seed_points,
                             const std::set<int>& feature_points,
                             const std::vector<Vector3>& bb_points,
                             RegularTriangulationNN& rt,
                             bool is_using_dilated_radius,
                             bool is_debug = false);
void generate_RT_dual_info(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,  // for winding filter inside/outside
    const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper,
    RegularTriangulationNN& rt, bool is_debug = false);

void update_dt(ThreeDimensionalShape& shape3D, bool is_debug);

class DtIO {
 public:
  // write valid ma spheres to file
  static void save_valid_medial_spheres(
      const std::vector<MVertex>& all_medial_spheres,
      const std::vector<int>& valid_medial_spheres);

  static void read_valid_medial_spheres(
      std::vector<MVertex>& all_medial_spheres,
      std::vector<int>& valid_medial_spheres,
      std::map<int, int>& all_to_valid_medial_spheres);

  // save random seeds during CVT
  static void save_seeds_random(const std::string& output_path,
                                const std::vector<double>& seed_points,
                                const std::vector<double>& seed_normals);

  static void load_seeds_random(const std::string& input_path,
                                std::vector<double>& seed_points,
                                std::vector<double>& seed_normals);

  // save seeds after CVT
  static void save_seeds_after_CVT(const std::string& output_path,
                                   const std::vector<double>& seed_points,
                                   const std::vector<double>& seed_normals,
                                   const std::set<int>& feature_points);

  static void load_seeds(const std::string& input_path,
                         std::vector<double>& seed_points,
                         std::vector<double>& seed_normals,
                         std::set<int>& feature_points);

  // save seeds after CVT for power crust
  static void save_seeds_to_pts(const std::string& output_name,
                                const Parameters& params,
                                const std::vector<double>& seed_points,
                                bool is_normalize);
  // same as save_seeds_to_pts(), just in .ma format
  static void save_seeds_to_ma(const std::string& output_name,
                               const Parameters& params,
                               const std::vector<double>& seed_points);

  static void save_all_mat_vertices(
      const std::string& output_path,
      const std::vector<double>& all_weighted_mat_pts_cgal,
      const std::vector<bool>& is_locking);

  static void load_all_mat_vertices(
      const std::string& input_path,
      std::vector<double>& all_weighted_mat_pts_cgal,
      std::vector<bool>& is_locking);

};  // DtIO

}  // namespace matfp