// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "matfp/Common.h"
#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/CommonTypes.h"

namespace matfp {

void load_sample_and_lock_sharp_edges(
    const MeshWrapper& sf_mesh_wrapper, const double& ideal_length,
    const double& s_edge_ideal_eps_rel,
    const GEO::NearestNeighborSearch_var& lfs_kd_tree,
    const std::vector<double>& lfs_min_values, const double r_sample,
    const std::set<int>& corners_set, std::vector<double>& seed_points,
    std::map<int, int>& seeds_map, std::set<int>& point_is_locked,
    std::vector<std::array<int, 2>>& s_edges,
    std::map<std::array<int, 2>, std::array<Vector3, 2>>& se_normals,
    std::vector<SpecialEdge>& sp_edges,
    std::map<std::array<int, 2>, int>& map_to_sp_edges, bool is_debug);

// will init tan_cc_lines
void load_concave_edges(const MeshWrapper& sf_mesh_wrapper,
                        std::vector<double>& seed_points,
                        std::map<int, int>& seeds_map,
                        std::vector<std::array<int, 2>>& cc_edges,
                        std::vector<SpecialEdge>& sp_edges,
                        std::vector<TangentConcaveLine>& tan_cc_lines,
                        std::map<std::array<int, 2>, int>& map_to_sp_edges, bool is_debug);

void load_and_lock_corners(const MeshWrapper& sf_mesh_wrapper,
                           std::vector<double>& seed_points,
                           std::map<int, int>& seeds_map,
                           std::set<int>& point_is_locked,
                           std::vector<int>& corners, bool is_debug);

void load_corner_special_edges(const std::vector<SpecialEdge>& sp_edges,
                               std::vector<int>& corners,
                               std::vector<SpecialCorner>& sp_corners, bool is_debug);

}  // namespace matfp