// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

void remove_medial_spheres_too_close_to_sharp_edges(
    const std::set<std::array<int, 2>>& se_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_using_dilated_radius,
    bool is_debug = false);

void prepare_all_medial_spheres(ThreeDimensionalShape& shape3D, bool is_debug);

bool add_or_delete_for_se(std::vector<MVertex>& all_medial_spheres,
                          std::set<std::array<int, 2>>& se_spheres,
                          std::map<int, int>& se_kd_tree_idx_to_se_tag,
                          GEO::NearestNeighborSearch_var& se_kd_tree,
                          std::vector<double>& se_kd_points,
                          bool is_check_updated_only,
                          bool is_using_dilated_radius, bool is_debug);

}  // namespace matfp