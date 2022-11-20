// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

// Shared functions
int get_num_adj_cc_parts(const MVertex& mat_A, const MVertex& mat_B);
int get_num_common_tan_normals(const MVertex& mat_A, const MVertex& mat_B);

// Main functions
void check_invalid_mat_edges(
    NonManifoldMesh& mat, const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres,  // not const on purpose
    std::set<std::array<int, 2>>& invalid_mat_edges, bool is_faster = true,
    bool is_debug = false);

void insert_internal_feature_spheres(
    const GEO::Mesh& sf_mesh,
    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
    const AABBWrapper& aabb_wrapper,
    const std::set<std::array<int, 2>>& invalid_mat_edges,
    const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_debug = false);
}  // namespace matfp
