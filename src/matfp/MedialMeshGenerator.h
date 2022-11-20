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
#include "matfp/ThreeDimensionalShape.h"
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/Types/RegularTriangulation.h"

namespace matfp {

using namespace std;
using namespace GEO;

void create_mat_from_RT(const std::vector<MVertex>& all_medial_spheres,
                        const std::vector<int>& valid_medial_spheres,
                        const RegularTriangulationNN& rt,
                        const std::string& mat_name, NonManifoldMesh& mat,
                        bool is_using_dilated_radius, bool is_debug);
}  // namespace matfp
