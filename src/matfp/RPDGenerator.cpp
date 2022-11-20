// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/RPDGenerator.h"

#include "ThreeDimensionalShape.h"
#include "matfp/Common.h"
#include "matfp/geogram/RPD.h"
#include <geogram/voronoi/RVD.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>

namespace matfp {

void update_rpd(ThreeDimensionalShape& shape3D, bool is_volumetric) {
  int dim = 4;
  igl::Timer timer;

  timer.start();
  // RPD
  if (!is_volumetric) {
    // TODO: rpd_seed_adj is not indicating correct connectivity
    //       only restricted cell that share the same edge should be connected
    matfp::generate_RPD_CGAL(shape3D.rt, shape3D.sf_mesh, shape3D.rpd_refined,
                             shape3D.rpd_merged, shape3D.rpd_seed_adj,
                             shape3D.rpd_vs_bisectors,
                             shape3D.rpd_cell_boundary_segments, is_volumetric);
  }
  // volumetric
  else {
    log_and_throw("ERROR: not handle volumetric now");
  }
  logger().info("RPD took {}s, is_volumetric {}", timer.getElapsedTimeInSec(),
                is_volumetric);
}

void export_RPD(const std::string rpd_name, const GEO::Mesh& rpd_mesh,
                const std::vector<int>& valid_medial_spheres,
                const std::vector<MVertex>& all_medial_spheres) {
  if (valid_medial_spheres.empty()) return;
  const GEO::Attribute<GEO::index_t> facet_region_attr(
      rpd_mesh.facets.attributes(), "region");

  std::string rpd_name_full =
      "../out/rpd/rpd_" + rpd_name + "_" + get_timestamp() + ".ma";
  logger().debug("start saving rpd to .ma file: {}", rpd_name_full);

  int nb_faces = 0;
  std::vector<int> faces;  // triangle faces
  std::vector<std::array<double, 3>> rgbs;
  for (int f = 0; f < rpd_mesh.facets.nb(); f++) {
    int seed_idx = facet_region_attr[f];
    if (seed_idx < 0 || seed_idx >= valid_medial_spheres.size()) continue;
    int f_nb_lv = rpd_mesh.facets.nb_vertices(f);
    if (f_nb_lv < 3) continue;
    int all_idx = valid_medial_spheres.at(seed_idx);
    const MVertex& mat_p = all_medial_spheres.at(all_idx);
    // store trianlge instead
    for (unsigned k = 1; k < f_nb_lv - 1; k++) {
      int v1 = rpd_mesh.facets.vertex(f, 0);
      int v2 = rpd_mesh.facets.vertex(f, k);
      int v3 = rpd_mesh.facets.vertex(f, k + 1);
      faces.push_back(v1);
      faces.push_back(v2);
      faces.push_back(v3);
      rgbs.push_back(mat_p.color_rgb);
      nb_faces++;
    }
  }
  if (faces.size() / 3 != nb_faces || rgbs.size() != nb_faces) {
    logger().error("nb_faces: {}, faces/3: {}, rgb: {}", nb_faces,
                   faces.size() / 3, rgbs.size());
    log_and_throw("ERROR!!!!");
  }

  std::ofstream fout(rpd_name_full);
  fout << rpd_mesh.vertices.nb() << " " << 0 << " " << nb_faces << std::endl;
  for (int v = 0; v < rpd_mesh.vertices.nb(); v++) {
    const GEO::vec3 p = rpd_mesh.vertices.point(v);
    fout << "v " << p[0] << " " << p[1] << " " << p[2] << " " << 0 << std::endl;
  }
  for (int f = 0; f < nb_faces; f++) {
    fout << "f ";
    fout << faces[f * 3] << " " << faces[f * 3 + 1] << " " << faces[f * 3 + 2]
         << " ";
    fout << rgbs[f][0] << " " << rgbs[f][1] << " " << rgbs[f][2] << std::endl;
  }
  fout.close();
}

// restricted power diagram
void generate_RPD_CGAL(
    RegularTriangulationNN& rt, GEO::Mesh& input, GEO::Mesh& rpd_refined,
    GEO::Mesh& rpd_merged,
    std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_seed_adj,
    std::map<GEO::index_t, std::set<GEO::index_t>>& rpd_vs_bisectors,
    std::map<GEO::index_t, std::set<std::array<GEO::index_t, 2>>>&
        rpd_cell_boundary_segments,
    const bool& volumetric) {
  if (volumetric) {
    log_and_throw("ERROR: not handle volumetric now");
  }

  logger().debug("Creating restricted power diagram...");
  if (rt.get_nb_vertices() < 1) {
    logger().error("[RPD] we have 0 vertices in RT, skip RPD");
    rpd_merged.clear();
    return;
  }

  matfp::RestrictedPowerDiagram* RPD =
      matfp::RestrictedPowerDiagram::create(&rt, &input);
  RPD->set_exact_predicates(true);
  RPD->set_volumetric(volumetric);

  rpd_refined.clear();
  rpd_seed_adj.clear();

  // compute in parallel
  // parallel enabled by default
  RPD->compute_RPD(rpd_refined, &rpd_seed_adj, &rpd_vs_bisectors, 3,
                   true, /*cell_borders_only*/
                   true, /*integration_simplices*/
                   false /*is_parallel*/
  );
  rpd_refined.assert_is_valid();
  print_mesh_info(rpd_refined);
}

}  // namespace matfp