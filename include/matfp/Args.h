// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <string>

namespace matfp {

// Global arguments controlling the behavior of MATFP
struct Args {
  // Initial target edge-length at every vertex (in % of the bbox diagonal)
  double initial_edge_len_rel = 1 / 20.0;

  // Input mesh path
  std::string input_surface_path = "";

  //////////////////////
  //  MATFP  options  //
  //////////////////////

  // When the surface is densely samples, we
  // use ds% of number of triangle centroids as our surface seeds,
  // easy for debugging since seeds are constants.
  double downsample_percentage = -1;

  // When downsample_percentage = -1, we sample random
  // surface seeds according to LFS. Seeds are randomly
  // sampled everytime so this is not easy for debugging.
  // Smaller the denser.
  double rsample = 0.9;

  // For adding medial spheres that tangent to concave lines using
  // ball-shrinking algo. We need to specify two parameters:
  // 1. cc_len_eps:
  //    defines the length between to pin points on concave lines, scaled in
  //    [0,10] (all models are normalized to [0, 10])
  // 2. cc_normal_eps:
  //    define the angle between two normals given a pin point, scaled in
  //    [0,360]
  // Given a pin point and a normal, we can add a sphere using ball-shrinking
  // algo. Concave lines can be better preserved when cc_len_eps is smaller
  // (more medial spheres inserted), but it also means more processing time.
  double cc_len_eps = 0.03;
  double cc_normal_eps = 10;

  //////////////////////
  // MATFP_PRE options//
  //////////////////////

  // Number of subdivison when input surface is too sparse
  int num_subdivide = 0;

  // Normalize it to [0, 10], please always do
  bool is_normalize = true;

  // Threshold for detecting sharp/concave edges and corners
  double thres_concave = 0.18;  // smaller more senstaive
  double thres_convex = 30.;    // smaller more senstaive

  // Flag for saving models as .geogram or not
  bool is_save_model = false;

  // Flag for saving scaled input model in .ply
  bool is_save_scaled_input = false;
};

}  // namespace matfp
