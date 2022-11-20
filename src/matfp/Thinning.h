// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <map>
#include <queue>
#include <unordered_set>
#include <vector>

#include "matfp/Logger.h"
#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

struct simple_pair {
  simple_pair(int _t, unsigned _idx0, unsigned _idx1, unsigned _tid) {
    type = _t;
    idx0 = _idx0;
    idx1 = _idx1;
    tid = _tid;
  }

  simple_pair(const simple_pair& _spr) {
    type = _spr.type;
    idx0 = _spr.idx0;
    idx1 = _spr.idx1;
    tid = _spr.tid;
  }

  simple_pair& operator=(const simple_pair& _spr) {
    type = _spr.type;
    idx0 = _spr.idx0;
    idx1 = _spr.idx1;

    return *this;
  }

  void print_info() const {
    logger().debug("----- simple pair, type {}, idx0 {}, idx1 {}, tid: {}",
                   type, idx0, idx1, tid);
  }

  // type = 2: tet-face pair, idx0 is tet id, idx1 is face id
  // type = 1: face-edge pair, idx0 is face id, idx1 is edge id
  // type = 0: edge-vert pair, idx0 is edge id, idx1 is vert id (no use)
  int type;
  unsigned idx0;
  unsigned idx1;

  // we only prune simple pairs from given tets, not other MAT elemetns
  // so here we store extra info about the tet for debugging
  unsigned tid;
};

typedef std::pair<double, int> face_importance;  // (importance, fid)

void prune(const std::vector<int>& valid_medial_spheres,
           const std::vector<MVertex>& all_medial_spheres, NonManifoldMesh& mat,
           double imp_thres = -1, bool is_sort_randomly = false,
           bool is_debug = false);
void load_all_mat_face_importance_globally(NonManifoldMesh& mat, bool is_debug);

}  // namespace matfp