// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/FeaturePreservation.h"

#include "AABBWrapper.h"
#include "IterateSpheres.h"
#include "matfp/InscribedSpheres.h"
#include "matfp/InternalFeatureAddition.h"

namespace matfp {

void search_external_feature(const std::vector<MVertex>& all_medial_spheres,
                             NonManifoldMesh& mat) {
  // save external medial feature
  mat.mat_extf_edges.clear();
  for (int v = 0; v < mat.numVertices; v++) {
    const auto& v_tmp = *(mat.vertices[v].second);
    const auto& mat_v_tmp = all_medial_spheres.at(v_tmp.all_idx);
    if (mat_v_tmp.is_on_s_edge) {
      for (const auto& adj_tag : mat_v_tmp.se_adj_se) {
        const auto& mat_v_adj = all_medial_spheres.at(adj_tag);
        if (mat_v_adj.valid_idx == -1) continue;
        std::array<int, 2> edge = {{v, mat_v_adj.valid_idx}};
        std::sort(edge.begin(), edge.end());
        mat.mat_extf_edges.insert(edge);
      }
    }
  }
  logger().debug("[Medial Feature] marked {} edges as mat external features",
                 mat.mat_extf_edges.size());
}
}  // namespace matfp