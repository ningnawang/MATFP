// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "matfp/ThreeDimensionalShape.h"

namespace matfp {

void generate_lfs(ThreeDimensionalShape& shape3D, bool is_debug = false);

}  // namespace matfp