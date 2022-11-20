// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "AABBWrapper.h"
#include "ThreeDimensionalShape.h"
#include "Triangulation.h"
#include "matfp/Types/DelaunayTriangulation.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

void inout_filter_raw(std::vector<MVertex>& all_medial_spheres,
                      const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
                      std::vector<double>& winding_num,
                      std::vector<bool>& inout, DelaunayTriangulation& dt);

// used during updating spheres to ideal positions
// in UpdateSpheres.cpp, some spheres might be pushed to
// outside of the object, so we will filter these out
void inout_filter_raw_simple(std::vector<MVertex>& all_medial_spheres,
                             const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F);
void inout_filter_vector_points(const std::vector<Vector3>& points,
                                const Eigen::MatrixXd& V,
                                const Eigen::MatrixXi& F,
                                std::vector<bool>& are_outside);

void convertSurface(const std::vector<Vector3>& input_vertices,
                    const std::vector<Vector3i>& input_faces,
                    Eigen::MatrixXd& V, Eigen::MatrixXi& F);

}  // namespace matfp