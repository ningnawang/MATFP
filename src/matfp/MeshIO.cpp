// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>
#include <igl/upsample.h>
#include <igl/write_triangle_mesh.h>
#include <matfp/Common.h>
#include <matfp/Logger.h>
#include <matfp/MeshIO.h>

namespace matfp {

// Here we assume mesh has been preprocessed
// with feature detected and stored in .geogram format (attributes)
bool MeshIO::load_mesh_from_geogram(const std::string& path, GEO::Mesh& input,
                                    std::vector<Vector3>& input_vertices,
                                    std::vector<Vector3i>& input_faces) {
  logger().debug("Loading mesh at {}...", path);
  if (getFileExtension(path) != "geogram") {
    log_and_throw("Please use mesh format as .geogram");
    return false;
  }
  input.clear(false, false);

  GEO::InputGeoFile geo_file(path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  // already called input.facets.connect() ???
  const bool ok = GEO::mesh_load(geo_file, input, flags);
  if (!ok) return false;

  MeshIO::init_vectors_from_mesh(input, input_vertices, input_faces);
  return ok;
}

// load points and faces from input
void MeshIO::init_vectors_from_mesh(const GEO::Mesh& input,
                                    std::vector<Vector3>& points,
                                    std::vector<Vector3i>& faces) {
  points.clear();
  points.resize(input.vertices.nb());
  for (size_t i = 0; i < points.size(); i++)
    points[i] << (input.vertices.point(i))[0], (input.vertices.point(i))[1],
        (input.vertices.point(i))[2];

  faces.clear();
  faces.resize(input.facets.nb());
  for (size_t i = 0; i < faces.size(); i++)
    faces[i] << input.facets.vertex(i, 0), input.facets.vertex(i, 1),
        input.facets.vertex(i, 2);
}

void MeshIO::export_sf_mesh_ma(const std::string sf_name,
                               const GEO::Mesh& sf_mesh) {
  std::string sf_name_full =
      "../out/input_scaled/input_" + sf_name + "_" + get_timestamp() + ".ma";
  logger().debug("start saving scaled input mesh to .ma file: {}",
                 sf_name_full);

  std::ofstream fout(sf_name_full);
  fout << sf_mesh.vertices.nb() << " " << 0 << " " << sf_mesh.facets.nb()
       << std::endl;
  for (int v = 0; v < sf_mesh.vertices.nb(); v++) {
    const GEO::vec3 p = sf_mesh.vertices.point(v);
    fout << "v " << p[0] << " " << p[1] << " " << p[2] << " "
         << 0 /*radius=0*/ << std::endl;
  }

  for (int f = 0; f < sf_mesh.facets.nb(); f++) {
    fout << "f ";
    int f_nb_lv = sf_mesh.facets.nb_vertices(f);
    for (int lv = 0; lv < f_nb_lv; lv++) {
      fout << sf_mesh.facets.vertex(f, lv) << " ";
    }
    fout << std::endl;
  }
  fout.close();
}

}  // namespace matfp