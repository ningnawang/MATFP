// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "pre_guiwindow.h"

#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <polyscope/volume_mesh.h>

#include "pre_meshIO.h"
#include "sharp_feature_detection.h"

namespace pre_matfp {
PreGuiWindow *PreGuiWindow::instance_ = nullptr;

PreGuiWindow::~PreGuiWindow() {
  delete m_shape;
  instance_ = nullptr;
}

PreGuiWindow::PreGuiWindow(Args &args) {
  if (instance_ != nullptr) {
    log_and_throw("ERROR: PreGuiWindow instance is not nullptr!!");
  }
  instance_ = this;
  std::string path = args.input_surface_path;
  std::string output_name = getFileName(args.input_surface_path, false);
  m_shape = new Shape3D();

  // Load mesh and detect features
  if (!load_mesh_and_preprocess(args, m_shape->sf_mesh)) {
    logger().error("Unable to load mesh at {}", path);
    return;
  }

  if (args.is_save_model) {
    save_mesh(output_name, m_shape->sf_mesh);
  }

  if (args.is_save_scaled_input) {
    save_scaled_input(output_name, m_shape->sf_mesh, false /*is_off*/);
  }
}

void PreGuiWindow::show() {
  ////////// Options
  // polyscope would make the tet in the center all the time
  // so we have no idea where the tet acutal location is
  // that's why we disabled polyscope::options::autocenterStructures in show()
  /* disble auto-center becaue we want to draw single tet for debugging */
  // polyscope::options::autocenterStructures = true;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;
  polyscope::view::style = polyscope::view::NavigateStyle::Free;
  polyscope::options::verbosity = 2;
  polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;

  // Initialize polyscope
  polyscope::init();

  // Register the mesh with Polyscope
  instance_->show_input_mesh(m_shape->sf_mesh);

  // Add the callback
  polyscope::state::userCallback = PreGuiWindow::callbacks;

  // Show the gui
  polyscope::show();
}

void PreGuiWindow::callbacks() {}

void PreGuiWindow::show_input_mesh(const GEO::Mesh &sf_mesh) {
  std::vector<Vector3> input_vertices;
  std::vector<Vector3i> input_faces;

  std::vector<Vector3> corner_vertices;
  std::set<std::array<int, 2>> s_edges, cc_edges;
  GEO::Attribute<int> attr_corners(sf_mesh.vertices.attributes(), "corner");
  GEO::Attribute<int> attr_se(sf_mesh.edges.attributes(), "se");
  GEO::Attribute<int> attr_cce(sf_mesh.edges.attributes(), "cce");

  for (auto v = 0; v < sf_mesh.vertices.nb(); v++) {
    input_vertices.push_back(to_eigen(sf_mesh.vertices.point(v)));
    if (attr_corners[v] == 1) {
      corner_vertices.push_back(input_vertices[v]);
    }
  }
  for (auto f = 0; f < sf_mesh.facets.nb(); f++) {
    Vector3i face;
    for (auto lv = 0; lv < 3; lv++) {
      face[lv] = sf_mesh.facets.vertex(f, lv);
    }
    input_faces.push_back(face);
  }
  for (int e = 0; e < sf_mesh.edges.nb(); e++) {
    std::array<int, 2> edge;
    for (int lv = 0; lv < 2; ++lv) {
      edge[lv] = sf_mesh.edges.vertex(e, lv);
    }
    if (attr_se[e] == 1) {
      s_edges.insert(edge);
    } else if (attr_cce[e] == 1) {
      cc_edges.insert(edge);
    }
  }

  // re-register
  auto input_mesh =
      polyscope::registerSurfaceMesh("input mesh", input_vertices, input_faces);
  input_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Cull);
  polyscope::registerPointCloud("Corners", corner_vertices);
  convert_to_show_edges(input_vertices, s_edges, "Sharp Edges");
  convert_to_show_edges(input_vertices, cc_edges, "Concave Edges");
}

void PreGuiWindow::convert_to_show_edges(
    const std::vector<Vector3> &old_pos,
    const std::set<std::array<int, 2>> &old_edges, std::string curve_name) {
  if (old_pos.empty() || old_edges.empty()) {
    // todo: remove existing edge
    return;
  }

  int new_idx = 0;
  std::map<int, int> old2new;
  std::vector<Vector3> new_pos;
  for (const auto &one_old_edge : old_edges) {
    for (const auto &v : one_old_edge) {
      if (old2new.find(v) != old2new.end()) {
        continue;
      }
      old2new[v] = new_idx;
      new_pos.push_back(old_pos[v]);
      new_idx++;
    }
  }

  std::set<std::array<int, 2>> new_edges;
  for (const auto &one_old_edge : old_edges) {
    std::array<int, 2> one_new_edge = {
        {old2new[one_old_edge[0]], old2new[one_old_edge[1]]}};
    std::sort(one_new_edge.begin(), one_new_edge.end());
    new_edges.insert(one_new_edge);
  }

  polyscope::registerCurveNetwork(curve_name, new_pos, new_edges);
}

}  // namespace pre_matfp