// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "pre_meshIO.h"

#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <igl/Timer.h>
#include <igl/remove_duplicate_vertices.h>
#include <igl/unique_rows.h>
#include <igl/upsample.h>
#include <igl/writeOFF.h>
#include <igl/writePLY.h>
#include <matfp/Logger.h>

// to use logger
using namespace matfp;

namespace pre_matfp {

std::string getFileName(std::string filePath, bool withExtension,
                        char seperator) {
  std::string filename_ext =
      filePath.substr(filePath.find_last_of(seperator) + 1);
  if (withExtension) return filename_ext;
  size_t lastindex = filename_ext.find_last_of(".");
  return filename_ext.substr(0, lastindex);
}

void save_mesh(const std::string sf_name, GEO::Mesh& mesh) {
  std::string output_path = "../out/mesh/mesh_" + sf_name + ".geogram";
  GEO::OutputGeoFile geo_file(output_path);
  GEO::MeshIOFlags flags;
  flags.set_element(GEO::MeshElementsFlags::MESH_ALL_ELEMENTS);
  flags.set_attributes(GEO::MeshAttributesFlags::MESH_ALL_ATTRIBUTES);
  if (!GEO::mesh_save(mesh, geo_file, flags)) {
    logger().error("Unable to save file at {}", output_path);
    return;
  }
  logger().debug("Done saving file {}", output_path);
}

// helper function for save_scaled_input()
void save_scaled_input(const std::string sf_name, const GEO::Mesh& sf_mesh,
                       bool is_off) {
  std::string sf_name_full = "../out/input_scaled/input_scaled_" + sf_name /*+
                             "_" + get_timestamp()*/
      ;

  // // save 0001510_partstudio_00_model_ste_00_2048.obj -> 1510.ply
  // std::string name_partial = sf_name.substr(0, 7);
  // if (name_partial[2] == '0') {
  //   name_partial = name_partial.substr(3);
  // } else {
  //   name_partial = name_partial.substr(2);
  // }
  // std::string sf_name_full = "../out/input_scaled/input_scaled_" +
  // name_partial;

  if (is_off) {
    sf_name_full += ".off";
  } else {
    sf_name_full += ".ply";
  }
  logger().debug("start saving scaled input mesh to .off file: {}",
                 sf_name_full);

  MatrixXs V(sf_mesh.vertices.nb(), 3);
  Eigen::MatrixXi F(sf_mesh.facets.nb(), 3);
  for (int i = 0; i < sf_mesh.vertices.nb(); i++) {
    V.row(i) = to_eigen(sf_mesh.vertices.point(i));
  }
  for (int i = 0; i < sf_mesh.facets.nb(); i++) {
    for (int j = 0; j < sf_mesh.facets.nb_vertices(i); j++)
      F(i, j) = sf_mesh.facets.vertex(i, j);
  }

  if (is_off)
    igl::writeOFF(sf_name_full, V, F);
  else
    igl::writePLY(sf_name_full, V, F);
}

// load points and faces from input
// will call input.facets.connect()
void init_vectors_from_mesh(const GEO::Mesh& input,
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

void reorder_mesh_from_vector(std::vector<Vector3>& points,
                              std::vector<Vector3i>& faces, GEO::Mesh& input) {
  logger().debug("Reordering mesh from internal data ...");
  input.clear(false, false);

  // Setup vertices
  input.vertices.create_vertices(points.size());
  for (int i = 0; i < input.vertices.nb(); ++i) {
    GEO::vec3& p = input.vertices.point(i);
    p[0] = points[i](0);
    p[1] = points[i](1);
    p[2] = points[i](2);
  }

  // Setup faces
  input.facets.create_triangles(faces.size());
  for (int c = 0; c < input.facets.nb(); ++c) {
    for (int lv = 0; lv < 3; ++lv) {
      input.facets.set_vertex(c, lv, faces[c](lv));
    }
  }

  // Setup edges
  std::vector<std::array<int, 2>> edges;
  for (int i = 0; i < faces.size(); i++) {
    const auto& f = faces[i];
    for (int j = 0; j < 3; j++) {
      std::array<int, 2> e = {{f[j], f[(j + 1) % 3]}};
      if (e[0] > e[1]) std::swap(e[0], e[1]);
      edges.push_back(e);
    }
  }
  vector_unique(edges);
  input.edges.create_edges(edges.size());
  for (int e = 0; e < edges.size(); e++) {
    for (int lv = 0; lv < 2; ++lv) {
      input.edges.set_vertex(e, lv, edges[e][lv]);
    }
  }

  GEO::mesh_reorder(input, GEO::MESH_ORDER_MORTON);
  // we did not setup the adjacent info till now
  input.facets.connect();

  // reload points and faces
  init_vectors_from_mesh(input, points, faces);
}

void normalize_mesh(const GEO::Mesh& mesh, std::vector<Vector3>& points,
                    std::vector<Vector3i>& faces) {
  const GEO::vec3& p0 = mesh.vertices.point(0);
  Vector3 min(p0[0], p0[1], p0[2]);
  Vector3 max(p0[0], p0[1], p0[2]);
  int size_max = 10;

  for (int i = 0; i < (int)mesh.vertices.nb(); ++i) {
    const GEO::vec3& p = mesh.vertices.point(i);
    for (int j = 0; j < 3; j++) {
      min(j) = std::min(min(j), p[j]);
      max(j) = std::max(max(j), p[j]);
    }
  }
  logger().debug("Rescaling & Centering");
  logger().debug("min corner:[{},{},{}]", min[0], min[1], min[2]);
  logger().debug("max corner:[{},{},{}]", max[0], max[1], max[2]);

  double size =
      std::max(max[0] - min[0], std::max(max[1] - min[1], max[2] - min[2])) /
      size_max;
  double xcenter = (max[0] + min[0]) * 0.5;
  double ycenter = (max[1] + min[1]) * 0.5;
  double zcenter = (max[2] + min[2]) * 0.5;

  points.clear();
  points.resize(mesh.vertices.nb());
  for (size_t i = 0; i < points.size(); i++)
    points[i] << ((mesh.vertices.point(i))[0] - xcenter) / size,
        ((mesh.vertices.point(i))[1] - ycenter) / size,
        ((mesh.vertices.point(i))[2] - zcenter) / size;

  faces.clear();
  faces.resize(mesh.facets.nb());
  for (size_t i = 0; i < faces.size(); i++)
    faces[i] << mesh.facets.vertex(i, 0), mesh.facets.vertex(i, 1),
        mesh.facets.vertex(i, 2);
}

bool remove_duplicates(std::vector<Vector3>& input_vertices,
                       std::vector<Vector3i>& input_faces) {
  MatrixXs V_tmp(input_vertices.size(), 3), V_in;
  Eigen::MatrixXi F_tmp(input_faces.size(), 3), F_in;
  for (int i = 0; i < input_vertices.size(); i++)
    V_tmp.row(i) = input_vertices[i];
  for (int i = 0; i < input_faces.size(); i++) F_tmp.row(i) = input_faces[i];

  //
  Eigen::VectorXi IV, _;
  igl::remove_duplicate_vertices(V_tmp, F_tmp, SCALAR_ZERO, V_in, IV, _, F_in);
  //
  for (int i = 0; i < F_in.rows(); i++) {
    int j_min = 0;
    for (int j = 1; j < 3; j++) {
      if (F_in(i, j) < F_in(i, j_min)) j_min = j;
    }
    if (j_min == 0) continue;
    int v0_id = F_in(i, j_min);
    int v1_id = F_in(i, (j_min + 1) % 3);
    int v2_id = F_in(i, (j_min + 2) % 3);
    F_in.row(i) << v0_id, v1_id, v2_id;
  }
  F_tmp.resize(0, 0);
  Eigen::VectorXi IF;
  igl::unique_rows(F_in, F_tmp, IF, _);
  F_in = F_tmp;
  //
  if (V_in.rows() == 0 || F_in.rows() == 0) return false;

  logger().info("remove duplicates: ");
  logger().info("#v: {} -> {}", input_vertices.size(), V_in.rows());
  logger().info("#f: {} -> {}", input_faces.size(), F_in.rows());

  input_vertices.resize(V_in.rows());
  input_faces.clear();
  input_faces.reserve(F_in.rows());
  for (int i = 0; i < V_in.rows(); i++) input_vertices[i] = V_in.row(i);
  for (int i = 0; i < F_in.rows(); i++) {
    if (F_in(i, 0) == F_in(i, 1) || F_in(i, 0) == F_in(i, 2) ||
        F_in(i, 2) == F_in(i, 1))
      continue;
    if (i > 0 && (F_in(i, 0) == F_in(i - 1, 0) &&
                  F_in(i, 1) == F_in(i - 1, 2) && F_in(i, 2) == F_in(i - 1, 1)))
      continue;
    // check area
    Vector3 u = V_in.row(F_in(i, 1)) - V_in.row(F_in(i, 0));
    Vector3 v = V_in.row(F_in(i, 2)) - V_in.row(F_in(i, 0));
    Vector3 area = u.cross(v);
    if (area.norm() / 2 <= SCALAR_ZERO) {
      logger().error(
          "ERROR: face {} area is too small, may impact our calculation", i);
      // continue;
    }
    input_faces.push_back(F_in.row(i));
  }

  return true;
}

bool remove_duplicates(MatrixXs& V_tmp, Eigen::MatrixXi& F_tmp,
                       std::vector<Vector3>& input_vertices,
                       std::vector<Vector3i>& input_faces) {
  MatrixXs V_in;
  Eigen::MatrixXi F_in;
  Eigen::VectorXi IV, _;
  igl::remove_duplicate_vertices(V_tmp, F_tmp, SCALAR_ZERO, V_in, IV, _, F_in);
  //
  for (int i = 0; i < F_in.rows(); i++) {
    int j_min = 0;
    for (int j = 1; j < 3; j++) {
      if (F_in(i, j) < F_in(i, j_min)) j_min = j;
    }
    if (j_min == 0) continue;
    int v0_id = F_in(i, j_min);
    int v1_id = F_in(i, (j_min + 1) % 3);
    int v2_id = F_in(i, (j_min + 2) % 3);
    F_in.row(i) << v0_id, v1_id, v2_id;
  }
  F_tmp.resize(0, 0);
  Eigen::VectorXi IF;
  igl::unique_rows(F_in, F_tmp, IF, _);
  F_in = F_tmp;
  if (V_in.rows() == 0 || F_in.rows() == 0) return false;

  logger().info("remove duplicates: ");
  logger().info("#v: {} -> {}", V_tmp.rows(), V_in.rows());
  logger().info("#f: {} -> {}", F_tmp.rows(), F_in.rows());

  input_vertices.resize(V_in.rows());
  input_faces.clear();
  input_faces.reserve(F_in.rows());
  for (int i = 0; i < V_in.rows(); i++) input_vertices[i] = V_in.row(i);
  for (int i = 0; i < F_in.rows(); i++) {
    if (F_in(i, 0) == F_in(i, 1) || F_in(i, 0) == F_in(i, 2) ||
        F_in(i, 2) == F_in(i, 1))
      continue;
    if (i > 0 && (F_in(i, 0) == F_in(i - 1, 0) &&
                  F_in(i, 1) == F_in(i - 1, 2) && F_in(i, 2) == F_in(i - 1, 1)))
      continue;
    // check area
    Vector3 u = V_in.row(F_in(i, 1)) - V_in.row(F_in(i, 0));
    Vector3 v = V_in.row(F_in(i, 2)) - V_in.row(F_in(i, 0));
    Vector3 area = u.cross(v);
    if (area.norm() / 2 <= SCALAR_ZERO) continue;
    input_faces.push_back(F_in.row(i));
  }

  return true;
}

void mesh_subdivision(std::vector<Vector3>& input_vertices,
                      std::vector<Vector3i>& input_faces,
                      const int num_subdivide) {
  logger().debug("start subdivisioning with num_subdivide: {}...",
                 num_subdivide);
  logger().debug("mesh before subdivision: #v: {}, #f: {}",
                 input_vertices.size(), input_faces.size());
  MatrixXs V(input_vertices.size(), 3), V_in;
  Eigen::MatrixXi F(input_faces.size(), 3), F_in;
  for (int i = 0; i < input_vertices.size(); i++) V.row(i) = input_vertices[i];
  for (int i = 0; i < input_faces.size(); i++) F.row(i) = input_faces[i];

  igl::upsample(V, F, V_in, F_in, num_subdivide);
  remove_duplicates(V_in, F_in, input_vertices, input_faces);

  logger().debug("mesh after subdivision: #v: {}, #f: {}",
                 input_vertices.size(), input_faces.size());
}
}  // namespace pre_matfp