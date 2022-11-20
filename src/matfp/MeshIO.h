// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include "matfp/Types/CommonTypes.h"

// /*
//  * Explicitly tell the moc NOT to parse openvdb header files.
//  * Details check README.md
// */
// #ifndef Q_MOC_RUN
// #undef Q_FOREACH
// #include <openvdb/openvdb.h>
// #endif

namespace matfp {

class MeshIO {
 public:
  static bool load_mesh_from_geogram(const std::string& path, GEO::Mesh& input,
                                     std::vector<Vector3>& input_vertices,
                                     std::vector<Vector3i>& input_faces);

  static void init_vectors_from_mesh(const GEO::Mesh& input,
                                     std::vector<Vector3>& points,
                                     std::vector<Vector3i>& faces);

  static void export_sf_mesh_ma(const std::string sf_name,
                                const GEO::Mesh& sf_mesh);
};

}  // namespace matfp
