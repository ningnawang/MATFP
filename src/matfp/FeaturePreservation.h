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
void search_external_feature(const std::vector<MVertex>& all_medial_spheres,
                             NonManifoldMesh& mat);

}  // namespace matfp