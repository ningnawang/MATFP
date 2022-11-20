// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "InscribedSpheres.h"
#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////////
// Push spheres to ideal positions
////////////////////////////////////////////////////////////////////////////////////

///////////////////////
// update mat spheres
// using RPD
/////////////////////////
void update_tan_elements_and_cc_parts_using_RDP(
    const GEO::Mesh& rpd_mesh,
    const std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_vs_bisectors,
    const std::vector<Vector3> ref_fs_normals, const bool is_volumetric,
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    const std::vector<int>& valid_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_update_tan_pls);

// for use
void update_spheres(const GEO::Mesh& sf_mesh,
                    const std::set<std::array<int, 2>>& ref_fs_pairs_not_cross,
                    const AABBWrapper& aabb_wrapper,
                    std::vector<MVertex>& all_medial_spheres, bool is_debug);

// remove low quality spheres
// 1. if contains only < 2 cc_part and no concave line, and normal variance is
// small
// 2. for non-feature sphere, all its cc_parts are around concave lines
// 3. extrude too much from any tangent plane
void clean_spheres(const std::vector<Vector3> ref_fs_normals,
                   std::vector<MVertex>& all_medial_spheres,
                   bool is_clean_extrude, bool is_clean_no_cc,
                   bool is_remove_dup);
}  // namespace matfp