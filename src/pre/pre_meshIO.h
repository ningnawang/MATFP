// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "matfp/Args.h"
#include "matfp/Logger.h"
#include "pre_types.h"

using namespace matfp;

namespace pre_matfp {
std::string getFileName(std::string filePath, bool withExtension = true,
                        char seperator = '/');
void save_mesh(const std::string sf_name, GEO::Mesh& mesh);
void save_scaled_input(const std::string sf_name, const GEO::Mesh& sf_mesh,
                       bool is_off);

// reorder and update
// will call input.facets.connect()
void reorder_mesh_from_vector(std::vector<Vector3>& points,
                              std::vector<Vector3i>& faces, GEO::Mesh& input);

void normalize_mesh(const GEO::Mesh& mesh, std::vector<Vector3>& points,
                    std::vector<Vector3i>& faces);

bool remove_duplicates(std::vector<Vector3>& input_vertices,
                       std::vector<Vector3i>& input_faces);

bool remove_duplicates(MatrixXs& V_tmp, Eigen::MatrixXi& F_tmp,
                       std::vector<Vector3>& input_vertices,
                       std::vector<Vector3i>& input_faces);

void mesh_subdivision(std::vector<Vector3>& input_vertices,
                      std::vector<Vector3i>& input_faces,
                      const int num_subdivide);

}  // namespace pre_matfp