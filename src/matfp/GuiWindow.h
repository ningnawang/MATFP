// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include "matfp/Args.h"
#include "matfp/ThreeDimensionalShape.h"

namespace matfp {

class GuiWindow {
 private:
  static GuiWindow* instance_;
  static constexpr double point_radius_rel = 0.0001;
  enum RestrictedType { INPUT, RVD_LFS, RVD, RPD, RPD_MERGED };
  enum MeshType { SURFACE, TET };

 private:
  ThreeDimensionalShape* m_shape3D = nullptr;
  Args* args = nullptr;

 public:
  GuiWindow(Args& args);
  ~GuiWindow();  // implemented in cpp

  // show polyscope window
  void show();

 private:
  static void callbacks();

  // callback functions
  void show_LFS(ThreeDimensionalShape* m_shape3D);
  void show_init_seeds(ThreeDimensionalShape* m_shape3D);
  void cal_rvd_or_rpd(ThreeDimensionalShape* m_shape3D,
                      const RestrictedType& type, const MeshType& mesh_type);

  // helper functions
  void show_restricted_diagram(ThreeDimensionalShape* m_shape3D,
                               const RestrictedType& type,
                               const MeshType& mesh_type);
  void show_restricted_mesh(const GEO::Mesh& restricted_mesh,
                            const RestrictedType& type,
                            const std::vector<int>& valid_medial_spheres,
                            const std::vector<MVertex>& all_medial_spheres);

  int given_mat_face_id = -1;
  void show_mat_simple(const NonManifoldMesh& mat, const RestrictedType& type,
                       int given_mat_face_id = -1,
                       bool is_show_unthin_trace = false);
  void convert_to_show_edges(const std::vector<Vector3>& old_pos,
                             const std::set<std::array<int, 2>>& old_edges,
                             std::string curve_name);

  int given_input_face_id = -1;
  void show_input_mesh(const GEO::Mesh& sf_mesh,
                       const int& given_input_face_id);

  double given_thinning_thres = 0.3;

 public:
  /////////////
  // Debug
  std::string mat_path;
  inline void set_mat_path(std::string path) { mat_path = path; }
};

}  // namespace matfp
