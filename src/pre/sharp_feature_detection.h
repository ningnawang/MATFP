// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include <unordered_set>
#include <vector>

#include "matfp/Args.h"
#include "matfp/Logger.h"
#include "pre_types.h"

using namespace matfp;

namespace pre_matfp {

bool load_mesh_and_preprocess(const Args& args, GEO::Mesh& input);

void find_feature_edges(const Args& args,
                        const std::vector<Vector3>& input_vertices,
                        const std::vector<Vector3i>& input_faces,
                        std::set<std::array<int, 2>>& s_edges,
                        std::set<std::array<int, 2>>& cc_edges,
                        std::set<int>& corners);

}  // namespace pre_matfp