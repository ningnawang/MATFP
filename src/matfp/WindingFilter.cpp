// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "WindingFilter.h"

#include <CGAL/Lazy_exact_nt.h>
#include <igl/winding_number.h>

#include <cmath>

#include "matfp/Logger.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/external/FastWindingNumber.hpp"

#define USE_FWN true

namespace matfp {

// detail checkout https://libigl.github.io/tutorial/#generalized-winding-number
// more inside -> 1, more outside -> 0
void inout_filter_raw(std::vector<MVertex>& all_medial_spheres,
                      const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                      std::vector<double>& winding_num,
                      std::vector<bool>& is_outside,
                      DelaunayTriangulation& dt) {
  // logger().debug("adding to C");
  int nb_pts = all_medial_spheres.size();
  Eigen::MatrixXd C(nb_pts, 3);  // all MAT
  for (int i = 0; i < nb_pts; i++) {
    Vector3& pos = all_medial_spheres[i].pos;
    for (int j = 0; j < 3; j++) {  // no need for weight info
      C(i, j) = pos[j];
    }
  }

  // logger().debug("winding_number");
  Eigen::VectorXd W;

#if USE_FWN
  matfp::fast_winding_number(V, F, C, W);
#else
  igl::winding_number(V, F, C, W);
#endif

  // logger().debug("output value");
  winding_num.clear();
  winding_num.resize(nb_pts);
  is_outside.clear();
  is_outside.resize(nb_pts);

  // determine inside/outside of cell
  int inside = 0;
  for (int i = 0; i < nb_pts; i++) {
    winding_num[i] = W(i);
    is_outside[i] = !(W(i) > 0.8);
    if (!is_outside[i]) inside++;

    // store info in all_medial_spheres
    all_medial_spheres[i].winding_num = W(i);
    all_medial_spheres[i].is_outside = !(W(i) > 0.8);
  }

  // store info in dt
  for (Finite_cells_iterator_dt fci = dt.finite_cells_begin();
       fci != dt.finite_cells_end(); fci++) {
    fci->info().winding_num = winding_num[fci->info().tag];
    fci->info().is_outside = is_outside[fci->info().tag];
  }

  logger().debug("filtered {}/{} inside mat vertices", inside, nb_pts);
  logger().info("In/out MAT filtering done!");
}

// used during updating spheres to ideal positions
// in UpdateSpheres.cpp, some spheres might be pushed to
// outside of the object, so we will filter these out
void inout_filter_raw_simple(std::vector<MVertex>& all_medial_spheres,
                             const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F) {
  int nb_spheres = all_medial_spheres.size();
  Eigen::MatrixXd C(nb_spheres, 3);  // all MAT
  for (int i = 0; i < nb_spheres; i++) {
    Vector3& pos = all_medial_spheres[i].pos;
    for (int j = 0; j < 3; j++) {  // no need for weight info
      C(i, j) = pos[j];
    }
  }

  // logger().debug("winding_number");
  Eigen::VectorXd W;

#if USE_FWN
  matfp::fast_winding_number(V, F, C, W);
#else
  igl::winding_number(V, F, C, W);
#endif

  // determine inside/outside of cell
  int inside = 0;
  for (int i = 0; i < nb_spheres; i++) {
    MVertex& mat_p = all_medial_spheres[i];

    if (mat_p.is_outside) continue;
    // do not update feature spheres
    if (mat_p.is_a_feature_sphere()) continue;

    // store winding info in all_medial_spheres
    mat_p.winding_num = W(i);
    mat_p.is_outside = !(W(i) > 0.9);
    if (!mat_p.is_outside) inside++;
  }

  logger().debug("[IN_OUT filter] keep {}/{} inside mat spheres", inside,
                 nb_spheres);
}

void inout_filter_vector_points(const std::vector<Vector3>& points,
                                const Eigen::MatrixXd& V,
                                const Eigen::MatrixXi& F,
                                std::vector<bool>& are_outside) {
  int nb_points = points.size();
  are_outside.clear();
  are_outside.resize(nb_points);
  Eigen::MatrixXd C(nb_points, 3);
  for (int i = 0; i < nb_points; i++) {
    const Vector3& pos = points[i];
    for (int j = 0; j < 3; j++) {
      C(i, j) = pos[j];
    }
  }

  Eigen::VectorXd W;
#if USE_FWN
  matfp::fast_winding_number(V, F, C, W);
#else
  igl::winding_number(V, F, C, W);
#endif

  // determine inside/outside of all points
  int inside = 0;
  for (int i = 0; i < nb_points; i++) {
    are_outside[i] = !(W(i) > 0.9);
    if (!are_outside[i]) inside++;
  }
  logger().debug("[IN_OUT filter] keep {}/{} inside points", inside, nb_points);
}

void convertSurface(const std::vector<Vector3>& input_vertices,
                    const std::vector<Vector3i>& input_faces,
                    Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  logger().debug("calling convertSurface ...");
  V.resize(input_vertices.size(), 3);
  Vector3 tmp;
  for (int i = 0; i < input_vertices.size(); i++) {
    Vector3 tmp = input_vertices[i];
    V(i, 0) = tmp[0];
    V(i, 1) = tmp[1];
    V(i, 2) = tmp[2];
  }

  F.resize(input_faces.size(), 3);
  for (int i = 0; i < input_faces.size(); i++) {
    for (int j = 0; j < 3; j++) F(i, j) = input_faces[i][j];
  }
  logger().debug("convertSurface done");
}

}  // namespace matfp
