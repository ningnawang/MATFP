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
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

////////////////////////////////////////////////////////////////////////////////////
// Common functions
// used by UpdateSpheres.h and IterateSpheres.h
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// DT related
////////////////////////////////////////////////////////////////////////////////////
enum TET_TYPE { T4_0 = 0, T3_1 = 1, T2_2 = 2, T2_1_1 = 3, T1_1_1_1 = 4 };

void collect_and_update_DT_feature_tets(
    const DelaunayTriangulation& dt, const std::vector<double>& seed_points,
    const std::vector<Vector3>& seed_normals_v2,
    std::vector<MVertex>& all_medial_spheres);

////////////////////////////////////////////////////////////////////////////////////
// NON-DT related
////////////////////////////////////////////////////////////////////////////////////
int remove_duplicated_medial_spheres(
    std::vector<MVertex>& sorted_partial_medial_spheres,
    std::vector<MVertex>& all_medial_spheres, bool is_override = false);

bool is_internal_medial_feature(const SphereType& type);

}  // namespace matfp
