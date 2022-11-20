// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

// To set the parameters related
#include <matfp/Types/CommonTypes.h>

#include <array>
#include <vector>

namespace matfp {

class Parameters {
 public:
  // initial target edge length at every vertex(in % of the box diagonal)
  double init_edge_length_rel = 1. / 20.;
  // epsilon presents the tolerence permited (in % of the box diagonal)
  double eps_rel = 1e-3;

  int init_nb_samples;
  double bbox_diag_length;
  double init_edge_length;
  double init_edge_length_2;
  double eps_input;

  double downsample_percentage = -1;

  // r-sample from Nina Amenta (check power crust)
  double r_sample = -1;

  // To make system robust
  bool is_using_dilated_radius = true;

  // Bbox
  bool is_bbox_stored = false;
  Vector3 bbox_min;
  Vector3 bbox_max;
  // for storing 8 bounding box vertices
  std::vector<Vector3> bb_points;

 public:
  bool init(double bbox_diag_l, double ds, double rs, double mesh_area = 0.) {
    bbox_diag_length = bbox_diag_l;
    downsample_percentage = ds;
    r_sample = rs;

    init_edge_length = bbox_diag_length * init_edge_length_rel;
    init_edge_length_2 = init_edge_length * init_edge_length;
    eps_input = bbox_diag_length * eps_rel;

    double ideal_area =
        std::sqrt(3) * 1. / 2. * init_edge_length * init_edge_length;
    init_nb_samples = std::ceil(mesh_area / ideal_area);  // estimation
    // std::cout << "init_nb_samples = " << init_nb_samples << std::endl;

    // std::cout << "bbox_diag_length = " << bbox_diag_length << std::endl;
    // std::cout << "init_edge_length = " << init_edge_length << std::endl;
    // std::cout << "downsample_percentage = " << downsample_percentage
    //           << std::endl;
    // std::cout << "r_sample = " << r_sample << std::endl;

    return true;
  }
};
}  // namespace matfp
