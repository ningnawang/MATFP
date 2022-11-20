// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"

namespace matfp {

void update_sphere_ss_param_using_seeds(
    const int nb_seeds, const std::set<int>& feature_points,
    std::vector<MVertex>& all_medial_spheres);

// will call shrink_sphere_wrapper()
void shrink_spheres(const GEO::Mesh& sf_mesh, const AABBWrapper& aabb_wrapper,
                    const std::vector<TangentConcaveLine>& tan_cc_lines,
                    std::vector<MVertex>& all_medial_spheres, bool is_debug);

// will call shrink_sphere()
bool shrink_sphere_wrapper(const GEO::Mesh& sf_mesh,
                           const AABBWrapper& aabb_wrapper, MVertex& mat_p,
                           bool is_setup_ss_param, bool is_check_cc,
                           bool is_debug);

// it performs better when cc_len_eps is smaller, and we do not need
// cc_normal_eps to be very small
void insert_spheres_for_concave_lines(
    const std::vector<TangentConcaveLine>& tan_cc_lines,
    std::vector<MVertex>& all_medial_spheres, double cc_len_eps = 2.,
    /*2 length, scaled in [0,10]*/
    double cc_normal_eps = 10. /*10 degree*/);

}  // namespace matfp