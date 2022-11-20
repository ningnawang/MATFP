// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/voronoi/RVD.h>
#include <igl/Timer.h>

#include <Eigen/Dense>

#include "matfp/Args.h"
#include "matfp/Common.h"
#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/CommonTypes.h"

namespace matfp {
using namespace GEO;

// Mesh processor
void mesh_remesh_split_sharp_edges(ThreeDimensionalShape& shape3D,
                                   bool is_debug);
//////////
void get_bb_corners(Parameters& params, const std::vector<double>& vertices,
                    Vector3& min, Vector3& max);

// sf_mesh is matching input_vertices & input_faces
void load_mesh_features(
    const GEO::Mesh& sf_mesh, const std::vector<Vector3>& input_vertices,
    const std::vector<Vector3i>& input_faces,
    std::vector<std::array<int, 2>>& s_edges,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals,
    std::set<std::array<int, 2>>& se_ref_fs_pairs,
    std::vector<std::array<int, 2>>& cc_edges, std::vector<int>& corners,
    std::map<int, std::unordered_set<int>>& conn_tris);

void init_input_faces_normals(const GEO::Mesh& input_mesh,
                              std::vector<Vector3>& input_fnormals);

}  // namespace matfp
