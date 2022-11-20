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

class PreGuiWindow {
 private:
  static PreGuiWindow* instance_;
  Shape3D* m_shape = nullptr;

 public:
  PreGuiWindow(Args& args);
  ~PreGuiWindow();  // implemented in cpp
  // show polyscope window
  void show();

  void show_input_mesh(const GEO::Mesh& sf_mesh);
  void convert_to_show_edges(const std::vector<Vector3>& old_pos,
                             const std::set<std::array<int, 2>>& old_edges,
                             std::string curve_name);

 private:
  static void callbacks();
};

}  // namespace pre_matfp