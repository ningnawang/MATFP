// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "matfp/AABBWrapper.h"
#include "matfp/Logger.h"
#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/Types/RegularTriangulation.h"

namespace matfp {

using namespace std;

// wrapper function
void update_rpd(ThreeDimensionalShape& shape3D,
                // bool is_refine_mat /*= false*/,
                // bool is_reload_all_weighted_mat_pts /*= true*/,
                // bool is_create_mat /*= true*/,
                // bool is_keep_previous_refine_result /*=false*/,
                bool is_volumetric = false);

// generate Restricted Power Diagram using CGAL
void generate_RPD_CGAL(
    RegularTriangulationNN& rt, GEO::Mesh& input, GEO::Mesh& rpd_refined,
    GEO::Mesh& rpd_merged,
    std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_seed_adj,
    std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_vs_bisectors,
    std::map<GEO::index_t, std::set<std::array<GEO::index_t, 2>>>&
        rpd_cell_boundary_segments,
    const bool& volumetric);

void export_RPD(const std::string rpd_name, const GEO::Mesh& rpd_mesh,
                const std::vector<int>& valid_medial_spheres,
                const std::vector<MVertex>& all_medial_spheres);

}  // namespace matfp
