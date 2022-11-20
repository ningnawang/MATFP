// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {
// Default parameters:
// alpha1 = 0.01;  // energy of distance to tangent point
// alpha2 = 1;     // energy of distance to tangent plane
// alpha3 = 1;     // energy of distance to concave line
bool iterate_sphere(const GEO::Mesh& sf_mesh,
                    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
                    const AABBWrapper& aabb_wrapper, MVertex& mat_p,
                    bool is_debug = false, double alpha1 = 0.01,
                    double alpha2 = 1, double alpha3 = 1,
                    const double break_threshold = SCALAR_ZERO_N4,
                    const int iteration_limit = 30);
}  // namespace matfp