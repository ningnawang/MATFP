// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "matfp/GuiWindow.h"

#include <igl/Timer.h>
#include <polyscope/curve_network.h>
#include <polyscope/point_cloud.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

#include "matfp/AABBWrapper.h"
#include "matfp/Common.h"
#include "matfp/FeaturePreservation.h"
#include "matfp/InternalFeatureAddition.h"
#include "matfp/IterateSpheres.h"
#include "matfp/LFS.h"
#include "matfp/Logger.h"
#include "matfp/MedialMeshGenerator.h"
#include "matfp/MedialSpheresProcessor.h"
#include "matfp/MeshIO.h"
#include "matfp/MeshProcessor.h"
#include "matfp/NonManifoldMesh/Nonmanifoldmesh.h"
#include "matfp/RPDGenerator.h"
#include "matfp/ShrinkSpheres.h"
#include "matfp/Thinning.h"
#include "matfp/Triangulation.h"
#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/OtherTypes.h"
#include "matfp/UpdateSpheres.h"
#include "matfp/WindingFilter.h"

namespace matfp {

GuiWindow *GuiWindow::instance_ = nullptr;

GuiWindow::~GuiWindow() {
  // we delete here since we allocate memory
  // (called new) for each of these variables
  delete m_shape3D;
  delete args;
  instance_ = nullptr;
}

GuiWindow::GuiWindow(Args &args) {
  if (instance_ != nullptr) {
    log_and_throw("ERROR: GuiWindow instance is not nullptr!!");
  }
  instance_ = this;
  std::string path = args.input_surface_path;
  this->args = &args;

  igl::Timer timer;
  bool skip_simplify = true;
  m_shape3D = new ThreeDimensionalShape;
  m_shape3D->mesh_name = getFileName(path, false);

  // Load mesh
  if (!MeshIO::load_mesh_from_geogram(path, m_shape3D->sf_mesh,
                                      m_shape3D->sf_mesh_wrapper.input_vertices,
                                      m_shape3D->sf_mesh_wrapper.input_faces)) {
    logger().error("Unable to load .geogram mesh at {}", path);
    return;
  }
  // Load sharp features
  matfp::load_mesh_features(
      m_shape3D->sf_mesh, m_shape3D->sf_mesh_wrapper.input_vertices,
      m_shape3D->sf_mesh_wrapper.input_faces, m_shape3D->s_edges,
      m_shape3D->se_normals, m_shape3D->se_ref_fs_pairs, m_shape3D->cc_edges,
      m_shape3D->corners, m_shape3D->sf_mesh_wrapper.conn_tris);
  matfp::init_input_faces_normals(m_shape3D->sf_mesh,
                                  m_shape3D->sf_mesh_wrapper.input_fnormals);
  // only load sharp edges, concave edges are not loaded
  // until load_concave_edges() is called
  m_shape3D->reload_ref_fs_pairs_not_cross();
  if (!m_shape3D->params.init(m_shape3D->get_sf_diag(),
                              args.downsample_percentage, args.rsample,
                              GEO::Geom::mesh_area(m_shape3D->sf_mesh))) {
    logger().error("Unable to init mesh parameters at {}", path);
    return;
  }

  // Init AABB
  // We do not reorder here because we assume MATFP_PRE has reordered,
  // see matfp_pre::load_mesh_and_preprocess()
  m_shape3D->aabb_wrapper.init_sf_mesh_and_tree(m_shape3D->sf_mesh,
                                                false /*is_reorder*/);
  m_shape3D->aabb_wrapper.init_feature_meshes_and_trees(
      m_shape3D->sf_mesh_wrapper.input_vertices, m_shape3D->s_edges,
      m_shape3D->cc_edges);

  // Fast winidng number
  matfp::convertSurface(m_shape3D->sf_mesh_wrapper.input_vertices,
                        m_shape3D->sf_mesh_wrapper.input_faces,
                        m_shape3D->sf_mesh_wrapper.VI,
                        m_shape3D->sf_mesh_wrapper.FI);
}

void GuiWindow::show() {
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
  instance_->show_input_mesh(m_shape3D->sf_mesh, -1 /*given_input_face_id*/);

  // Add the callback
  polyscope::state::userCallback = GuiWindow::callbacks;

  // Show the gui
  polyscope::show();
}

void GuiWindow::callbacks() {
  igl::Timer timer;
  timer.start();
  ImGui::PushItemWidth(100);

  if (ImGui::SmallButton("Stage1")) {
    matfp::generate_lfs(*(instance_->m_shape3D), false /*is_debug*/);
    // instance_->show_LFS(instance_->m_shape3D);
    matfp::mesh_remesh_split_sharp_edges(*(instance_->m_shape3D),
                                         false /*is_debug*/);
    instance_->show_init_seeds(instance_->m_shape3D);
    matfp::update_dt(*(instance_->m_shape3D), false /*is_debug*/);
    matfp::prepare_all_medial_spheres(*(instance_->m_shape3D),
                                      false /*is_debug*/);

    // it performs better when cc_len_eps is smaller, and we do not need
    // cc_normal_eps to be very small
    matfp::insert_spheres_for_concave_lines(
        instance_->m_shape3D->tan_cc_lines,
        instance_->m_shape3D->all_medial_spheres, instance_->args->cc_len_eps,
        instance_->args->cc_normal_eps);
    logger().info("Stage1 took {}s", timer.getElapsedTimeInSec());
  }

  ImGui::SameLine();
  if (ImGui::SmallButton("Stage2")) {
    matfp::shrink_spheres(
        instance_->m_shape3D->sf_mesh, instance_->m_shape3D->aabb_wrapper,
        instance_->m_shape3D->tan_cc_lines,
        instance_->m_shape3D->all_medial_spheres, false /*is_debug*/);
    matfp::update_spheres(instance_->m_shape3D->sf_mesh,
                          instance_->m_shape3D->ref_fs_pairs_not_cross,
                          instance_->m_shape3D->aabb_wrapper,
                          instance_->m_shape3D->all_medial_spheres,
                          false /*is_debug*/);
    logger().info("Stage2 took {}s", timer.getElapsedTimeInSec());
  }

  ImGui::SameLine();
  if (ImGui::SmallButton("Stage3")) {
    matfp::add_or_delete_for_se(
        instance_->m_shape3D->all_medial_spheres,
        instance_->m_shape3D->se_spheres,
        instance_->m_shape3D->se_kd_tree_idx_to_se_tag,
        instance_->m_shape3D->se_kd_tree, instance_->m_shape3D->se_kd_points,
        true /*is_check_updated_only*/,
        instance_->m_shape3D->params.is_using_dilated_radius,
        false /*is_debug*/);
    matfp::generate_RT_for_dual_PD(
        instance_->m_shape3D->all_medial_spheres,
        instance_->m_shape3D->valid_medial_spheres,
        instance_->m_shape3D->seed_points, instance_->m_shape3D->feature_points,
        instance_->m_shape3D->params.bb_points, instance_->m_shape3D->rt,
        instance_->m_shape3D->params.is_using_dilated_radius,
        false /*is_debug*/);
    matfp::generate_RT_dual_info(instance_->m_shape3D->sf_mesh_wrapper.VI,
                                 instance_->m_shape3D->sf_mesh_wrapper.FI,
                                 instance_->m_shape3D->sf_mesh,
                                 instance_->m_shape3D->aabb_wrapper,
                                 instance_->m_shape3D->rt, false /*is_debug*/);
    matfp::create_mat_from_RT(
        instance_->m_shape3D->all_medial_spheres,
        instance_->m_shape3D->valid_medial_spheres, instance_->m_shape3D->rt,
        instance_->m_shape3D->mesh_name, instance_->m_shape3D->mat_refined,
        instance_->m_shape3D->params.is_using_dilated_radius,
        false /*is_debug*/);
    instance_->show_mat_simple(instance_->m_shape3D->mat_refined,
                               RestrictedType::RPD);
    logger().info("Stage3 took {}s", timer.getElapsedTimeInSec());
  }

  ImGui::SameLine();
  if (ImGui::SmallButton("Stage4")) {
    for (int i = 0; i < 2; i++) {
      matfp::check_invalid_mat_edges(
          instance_->m_shape3D->mat_refined,
          instance_->m_shape3D->valid_medial_spheres,
          instance_->m_shape3D->all_medial_spheres,
          instance_->m_shape3D->mat_refined.invalid_mat_edges,
          true /*is_faster*/, false /*is_debug*/);
      // logger().info("Check_invalid took {}s", timer.getElapsedTimeInSec());
      matfp::insert_internal_feature_spheres(
          instance_->m_shape3D->sf_mesh,
          instance_->m_shape3D->ref_fs_pairs_not_cross,
          instance_->m_shape3D->aabb_wrapper,
          instance_->m_shape3D->mat_refined.invalid_mat_edges,
          instance_->m_shape3D->valid_medial_spheres,
          instance_->m_shape3D->all_medial_spheres, false /*is_debug*/);
      matfp::add_or_delete_for_se(
          instance_->m_shape3D->all_medial_spheres,
          instance_->m_shape3D->se_spheres,
          instance_->m_shape3D->se_kd_tree_idx_to_se_tag,
          instance_->m_shape3D->se_kd_tree, instance_->m_shape3D->se_kd_points,
          true /*is_check_updated_only*/,
          instance_->m_shape3D->params.is_using_dilated_radius,
          false /*is_debug*/);
      matfp::generate_RT_for_dual_PD(
          instance_->m_shape3D->all_medial_spheres,
          instance_->m_shape3D->valid_medial_spheres,
          instance_->m_shape3D->seed_points,
          instance_->m_shape3D->feature_points,
          instance_->m_shape3D->params.bb_points, instance_->m_shape3D->rt,
          instance_->m_shape3D->params.is_using_dilated_radius,
          false /*is_debug*/);
      matfp::generate_RT_dual_info(
          instance_->m_shape3D->sf_mesh_wrapper.VI,
          instance_->m_shape3D->sf_mesh_wrapper.FI,
          instance_->m_shape3D->sf_mesh, instance_->m_shape3D->aabb_wrapper,
          instance_->m_shape3D->rt, false /*is_debug*/);
      matfp::create_mat_from_RT(
          instance_->m_shape3D->all_medial_spheres,
          instance_->m_shape3D->valid_medial_spheres, instance_->m_shape3D->rt,
          instance_->m_shape3D->mesh_name, instance_->m_shape3D->mat_refined,
          instance_->m_shape3D->params.is_using_dilated_radius,
          false /*is_debug*/);
    }
    instance_->show_mat_simple(instance_->m_shape3D->mat_refined,
                               RestrictedType::RPD);
    logger().info("Stage4 took {}s", timer.getElapsedTimeInSec());
  }

  // Prune / Thinning
  // threshold: more faces will be prunes if bigger
  ImGui::InputDouble("threshold", &(instance_->given_thinning_thres));
  ImGui::SameLine();
  if (ImGui::SmallButton("Thinning")) {
    matfp::prune(instance_->m_shape3D->valid_medial_spheres,
                 instance_->m_shape3D->all_medial_spheres,
                 instance_->m_shape3D->mat_refined,
                 instance_->given_thinning_thres, false /*is_sort_randomly*/,
                 false /*is_debug*/);
    instance_->show_mat_simple(instance_->m_shape3D->mat_refined,
                               RestrictedType::RPD);
    logger().info("Thinning took {}s", timer.getElapsedTimeInSec());
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("show importance")) {
    matfp::load_all_mat_face_importance_globally(
        instance_->m_shape3D->mat_refined, false /*is_debug*/);
    instance_->show_mat_simple(instance_->m_shape3D->mat_refined,
                               RestrictedType::RPD);
  }

  // RPD
  if (ImGui::SmallButton("Calculate RPD")) {
    matfp::generate_RT_for_dual_PD(
        instance_->m_shape3D->all_medial_spheres,
        instance_->m_shape3D->valid_medial_spheres,
        instance_->m_shape3D->seed_points, instance_->m_shape3D->feature_points,
        instance_->m_shape3D->params.bb_points, instance_->m_shape3D->rt,
        instance_->m_shape3D->params.is_using_dilated_radius,
        false /*is_debug*/);
    instance_->cal_rvd_or_rpd(instance_->m_shape3D, RestrictedType::RPD,
                              MeshType::SURFACE);
    instance_->show_restricted_diagram(instance_->m_shape3D,
                                       RestrictedType::RPD, MeshType::SURFACE);
    // matfp::update_tan_elements_and_cc_parts_using_RDP(
    //     instance_->m_shape3D->rpd_refined,
    //     instance_->m_shape3D->rpd_vs_bisectors,
    //     instance_->m_shape3D->sf_mesh_wrapper.input_fnormals,
    //     instance_->m_shape3D->is_volumetric ? true : false, /*is_volumetric*/
    //     instance_->m_shape3D->tan_cc_lines,
    //     instance_->m_shape3D->valid_medial_spheres,
    //     instance_->m_shape3D->all_medial_spheres, true
    //     /*is_update_tan_pls*/);
    logger().info("Init_RPD took {}s", timer.getElapsedTimeInSec());
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("Save RPD")) {
    matfp::export_RPD(instance_->m_shape3D->mesh_name,
                      instance_->m_shape3D->rpd_refined,
                      instance_->m_shape3D->valid_medial_spheres,
                      instance_->m_shape3D->all_medial_spheres);
  }

  if (ImGui::SmallButton("save MAT .ply")) {
    MatIO::write_nmm_ply(instance_->m_shape3D->mat_refined.mat_name,
                         instance_->m_shape3D->mat_refined);
    logger().info("Saved MAT to {}.ply",
                  instance_->m_shape3D->mat_refined.mat_name);
  }
  ImGui::SameLine();
  if (ImGui::SmallButton("save MAT .ma")) {
    MatIO::export_nmm(instance_->m_shape3D->mat_refined.mat_name,
                      instance_->m_shape3D->mat_refined);
    logger().info("Saved MAT to {}.ma",
                  instance_->m_shape3D->mat_refined.mat_name);
  }

  ImGui::PopItemWidth();
}

void GuiWindow::show_input_mesh(const GEO::Mesh &sf_mesh,
                                const int &given_input_face_id) {
  std::vector<Vector3> input_vertices;
  std::vector<Vector3i> input_faces;

  std::vector<Vector3> corner_vertices;
  std::set<std::array<int, 2>> s_edges, cc_edges;
  GEO::Attribute<int> attr_corners(sf_mesh.vertices.attributes(), "corner");
  GEO::Attribute<int> attr_se(sf_mesh.edges.attributes(), "se");
  GEO::Attribute<int> attr_cce(sf_mesh.edges.attributes(), "cce");

  std::vector<bool> is_selected(sf_mesh.facets.nb(), false);
  for (auto v = 0; v < sf_mesh.vertices.nb(); v++) {
    input_vertices.push_back(to_eigen(sf_mesh.vertices.point(v)));
    if (attr_corners[v] == 1) {
      corner_vertices.push_back(input_vertices[v]);
    }
  }
  for (auto f = 0; f < sf_mesh.facets.nb(); f++) {
    if (given_input_face_id == f) {
      is_selected[f] = true;
    }
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
  if (given_input_face_id == -1) {
    auto input_mesh = polyscope::registerSurfaceMesh(
        "Input mesh", input_vertices, input_faces);
    input_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Cull);
    polyscope::registerPointCloud("Corners", corner_vertices);
    convert_to_show_edges(input_vertices, s_edges, "Sharp Edges");
    convert_to_show_edges(input_vertices, cc_edges, "Concave Edges");
  } else if (polyscope::hasSurfaceMesh("Input mesh")) {
    auto input_mesh = polyscope::getSurfaceMesh("Input mesh");
    input_mesh->addFaceScalarQuantity("is_selected", is_selected)
        ->setEnabled(true);
  }
}

void GuiWindow::show_LFS(ThreeDimensionalShape *m_shape3D) {
  const GEO::Mesh &sf_mesh = m_shape3D->sf_mesh;
  const GEO::NearestNeighborSearch_var lfs_kd_tree = m_shape3D->lfs_kd_tree;
  const double *lfs_min_values = m_shape3D->lfs_min_values.data();

  VectorXs lfs(sf_mesh.vertices.nb());
  for (int v = 0; v < sf_mesh.vertices.nb(); v++) {
    const double *ptr = sf_mesh.vertices.point_ptr(v);
    int lfs_idx = lfs_kd_tree->get_nearest_neighbor(ptr);
    double lfs_value = *(lfs_min_values + lfs_idx);
    lfs[v] = lfs_value;
  }
  polyscope::getSurfaceMesh("Input mesh")
      ->addVertexScalarQuantity("lfs", lfs, polyscope::DataType::MAGNITUDE);
}

void GuiWindow::show_init_seeds(ThreeDimensionalShape *m_shape3D) {
  int init_seeds_size = m_shape3D->init_seed_points.size() / 3;
  std::vector<Vector3> init_seeds(init_seeds_size);
  for (int i = 0; i < init_seeds_size; i++) {
    for (int j = 0; j < 3; j++) {
      init_seeds[i][j] = m_shape3D->init_seed_points[i * 3 + j];
    }
  }
  polyscope::registerPointCloud("Init seeds", init_seeds);
  polyscope::getPointCloud("Init seeds")
      ->setPointRadius(GuiWindow::point_radius_rel)
      ->setEnabled(false);

  // others
  // update face density: fd_rho and lfs_rho
  GEO::Attribute<double> facet_density_attr(
      m_shape3D->sf_mesh.facets.attributes(), "density");
  VectorXs facet_densitis(facet_density_attr.size());
  for (size_t i = 0; i < facet_density_attr.size(); i++) {
    facet_densitis[i] = facet_density_attr[i];
  }
  polyscope::getSurfaceMesh("Input mesh")
      ->addFaceScalarQuantity("face_density", facet_densitis,
                              polyscope::DataType::MAGNITUDE);
}

void GuiWindow::cal_rvd_or_rpd(ThreeDimensionalShape *m_shape3D,
                               const RestrictedType &type,
                               const MeshType &mesh_type) {
  if (type == RestrictedType::RVD) {
  } else if (type == RestrictedType::RPD) {
    if (mesh_type == MeshType::SURFACE)
      matfp::update_rpd(*m_shape3D);
    else if (mesh_type == MeshType::TET) {
      matfp::update_rpd(*m_shape3D, true);
    }
  }
}

void GuiWindow::show_restricted_diagram(ThreeDimensionalShape *m_shape3D,
                                        const RestrictedType &type,
                                        const MeshType &mesh_type) {
  // show surfaic RPD
  instance_->show_restricted_mesh(m_shape3D->rpd_refined, RestrictedType::RPD,
                                  m_shape3D->valid_medial_spheres,
                                  m_shape3D->all_medial_spheres);
  // others
  if (polyscope::hasSurfaceMesh("sample points"))
    polyscope::getPointCloud("sample points")->setEnabled(false);
  if (polyscope::hasSurfaceMesh("Input mesh"))
    polyscope::getSurfaceMesh("Input mesh")->setEnabled(false);
}

void GuiWindow::show_restricted_mesh(
    const GEO::Mesh &restricted_mesh, const RestrictedType &type,
    const std::vector<int> &valid_medial_spheres,
    const std::vector<MVertex> &all_medial_spheres) {
  const GEO::Attribute<GEO::index_t> facet_region_attr(
      restricted_mesh.facets.attributes(), "region");
  const GEO::Attribute<GEO::index_t> facet_ref_facet_attr(
      restricted_mesh.facets.attributes(), "ref_facet");
  // const GEO::Attribute<GEO::index_t> facet_ref_segment_attr(
  // 	restricted_mesh.facets.attributes(), "ref_segment"
  // );
  std::vector<int> all_mat_tags(restricted_mesh.facets.nb(), -1);
  std::vector<int> valid_mat_tags(restricted_mesh.facets.nb(), -1);
  std::vector<int> ref_facets(restricted_mesh.facets.nb(),
                              -1);  // original facet indices
  std::vector<int> ref_segments(restricted_mesh.facets.nb(),
                                -1);  // original facet indices

  std::vector<Vector3> r_mesh_vertices(restricted_mesh.vertices.nb());
  std::vector<VectorXi> r_mesh_faces(
      restricted_mesh.facets.nb());  // vs size might be > 3
  std::vector<std::array<double, 3>> f_color(restricted_mesh.facets.nb());
  std::map<int, std::array<double, 3>> seed_rgb;

  for (int v = 0; v < restricted_mesh.vertices.nb(); v++) {
    const GEO::vec3 p = restricted_mesh.vertices.point(v);
    r_mesh_vertices[v] << p[0], p[1], p[2];
  }

  double r = 0., g = 0., b = 0.;
  for (int f = 0; f < restricted_mesh.facets.nb(); f++) {
    int seed_idx = facet_region_attr[f];
    valid_mat_tags[f] = seed_idx;
    if (valid_medial_spheres.size() > 0 && seed_idx >= 0 &&
        seed_idx < valid_medial_spheres.size()) {
      all_mat_tags[f] = valid_medial_spheres.at(seed_idx);
      const MVertex &mat_p = all_medial_spheres.at(all_mat_tags[f]);
      f_color[f] = mat_p.color_rgb;
    } else {
      // if (seed_rgb.find(seed_idx) != seed_rgb.end())
      // {
      //     std::array<double, 3> &rgb = seed_rgb[facet_region_attr[f]];
      //     f_color[f] = {{rgb[0], rgb[1], rgb[2]}};
      // }
      // else
      // {
      //     matfp::color_addon_helper(r, g, b);
      //     std::array<double, 3> rgb = {{r, g, b}};
      //     seed_rgb.insert(
      //         std::pair<int, std::array<double, 3>>(seed_idx, rgb));
      //     f_color[f] = {{rgb[0], rgb[1], rgb[2]}};
      // }
    }

    int f_nb_v = restricted_mesh.facets.nb_vertices(f);
    if (f_nb_v < 3) {
      logger().debug("facet {} has < 3 vertices, mat valid_idx {}", f,
                     seed_idx);
      for (int lv = 0; lv < f_nb_v; lv++) {
        int v = restricted_mesh.facets.vertex(f, lv);
        logger().debug("face {} contain vertex {}", f, v);
      }
    }

    r_mesh_faces[f].resize(f_nb_v);
    for (int lv = 0; lv < f_nb_v; lv++) {
      int v = restricted_mesh.facets.vertex(f, lv);
      r_mesh_faces[f][lv] = v;
    }

    ref_facets[f] = facet_ref_facet_attr[f];
    // ref_segments[f] = facet_ref_segment_attr[f];
  }

  polyscope::SurfaceMesh *poly_mesh = nullptr;
  if (type == RestrictedType::RVD_LFS) {
    poly_mesh = polyscope::registerSurfaceMesh("RVD for LFS", r_mesh_vertices,
                                               r_mesh_faces);
  } else if (type == RestrictedType::RVD) {
    poly_mesh = polyscope::registerSurfaceMesh("RVD after CVT", r_mesh_vertices,
                                               r_mesh_faces);
  } else if (type == RestrictedType::RPD) {
    poly_mesh =
        polyscope::registerSurfaceMesh("RPD", r_mesh_vertices, r_mesh_faces);
  }
  if (poly_mesh != nullptr) {
    poly_mesh->addFaceColorQuantity("f_color", f_color)->setEnabled(true);
    poly_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Cull);
  }

  // store original facet indices
  poly_mesh->addFaceScalarQuantity("ref_facet:", ref_facets);
  // poly_mesh->addFaceScalarQuantity("ref_segment:", ref_segments);

  // store mat_p indices
  poly_mesh->addFaceScalarQuantity("mat_p valid_idx", valid_mat_tags);
  poly_mesh->addFaceScalarQuantity("mat_p all_idx", all_mat_tags);
}

void GuiWindow::show_mat_simple(const NonManifoldMesh &mat,
                                const RestrictedType &type,
                                int given_mat_face_id,
                                bool is_show_unthin_trace) {
  std::vector<Vector3> mat_pos(mat.numVertices);
  std::vector<double> mat_radius(mat.numVertices);
  std::vector<Vector3i> mat_faces(mat.numFaces, Vector3i(0, 0, 0));
  std::vector<bool> mat_faces_unthin(mat.numFaces, false);
  std::vector<bool> is_given_mat_face(
      mat.numFaces, false);  // set true when given_mat_face_id != -1

  // for showing dual segment of a RT face
  std::vector<bool> is_dual_seg_intersect(mat.numFaces, true);
  std::vector<bool> mat_faces_deleted(mat.numFaces, false);
  std::vector<Vector3> dual_segment_vs;
  std::vector<std::array<int, 2>> dual_segment_edge;
  // for showing intersection segment of dual_segment and input surface
  std::vector<Vector3> dual_intersect_segement_vs;
  std::vector<std::array<int, 2>> dual_intersect_segement_edge;
  std::vector<double> mat_faces_importance(mat.numFaces, 0.);

  // store mat vertices
  for (int i = 0; i < mat.numVertices; i++) {
    mat_pos[i] = mat.vertices[i].second->pos;
    mat_radius[i] = mat.vertices[i].second->radius;
  }

  // fetch all mat faces ids to show up
  std::vector<int> all_mat_faces_idx;
  if (is_show_unthin_trace) {
    all_mat_faces_idx = mat.unthin_trace;
  } else {
    for (int f = 0; f < mat.numFaces; f++) all_mat_faces_idx.push_back(f);
  }

  // store mat faces
  // need to enbale culling for drawing ma faces twice
  for (const int f : all_mat_faces_idx) {
    // if face is deleted, then keep 0 as vertex indices
    if (mat.faces[f].second->is_deleted) {
      continue;
    }
    mat_faces_importance[f] = mat.faces[f].second->importance;
    is_dual_seg_intersect[f] = mat.faces[f].second->is_dual_seg_intersect;
    mat_faces_deleted[f] = mat.faces[f].second->is_deleted;

    // mark if mat face is problematic (unthin that makes euler change)
    if (mat.unthin_faces.find(f) != mat.unthin_faces.end()) {
      mat_faces_unthin[f] = true;
      // mat_faces_unthin[mat.numFaces + f] = true;
    }

    // draw facets (counter-clockwise)
    int idx = 0;
    for (auto si = mat.faces[f].second->vertices_.begin();
         si != mat.faces[f].second->vertices_.end(); si++) {
      mat_faces[f][idx] = *si;
      idx++;
      Vector3 &pos = mat.vertices[*si].second->pos;
    }

    // // draw facets reversely (clockwise)
    // idx = 0;
    // for (auto si = mat.faces[f].second->vertices_.rbegin();
    //      si != mat.faces[f].second->vertices_.rend(); si++)
    // {
    //     mat_faces[mat.numFaces + f][idx] = *si;
    //     idx++;
    // }
  }

  // if this mat face is what we are looking for
  if (given_mat_face_id > 0 && given_mat_face_id < mat.numFaces) {
    int f = given_mat_face_id;
    const auto &matf = *(mat.faces[f].second);
    is_given_mat_face[f] = true;
    dual_segment_vs.push_back(matf.dual_segment[0]);
    dual_segment_vs.push_back(matf.dual_segment[1]);
    dual_segment_edge.push_back({{0, 1}});

    dual_intersect_segement_vs.push_back(matf.dual_intersections[0]);
    dual_intersect_segement_vs.push_back(matf.dual_intersections[1]);
    dual_intersect_segement_edge.push_back({{0, 1}});

    // print face info
    matf.print_mat_face();
    std::vector<int> f_vs_tags;
    for (const auto &vid : matf.vertices_) {
      f_vs_tags.push_back(mat.vertices[vid].second->all_idx);
    }
    logger().debug("mat face {} has vertices tags: {}", f, f_vs_tags);
  }

  // Register medial mesh
  polyscope::SurfaceMesh *medial_mesh = nullptr;
  if (is_show_unthin_trace) {
    medial_mesh =
        polyscope::registerSurfaceMesh("Medial mesh trace", mat_pos, mat_faces);
  } else if (type == RestrictedType::RVD) {
    medial_mesh =
        polyscope::registerSurfaceMesh("Medial mesh RVD", mat_pos, mat_faces);
  } else if (type == RestrictedType::RPD) {
    medial_mesh =
        polyscope::registerSurfaceMesh("Medial mesh RPD", mat_pos, mat_faces);
  }
  if (medial_mesh != nullptr) {
    medial_mesh->addFaceScalarQuantity("unthin", mat_faces_unthin);
    if (given_mat_face_id != -1) {
      medial_mesh->addFaceScalarQuantity("is_given_mat_face", is_given_mat_face)
          ->setEnabled(true);
      polyscope::registerCurveNetwork("MM dual segment", dual_segment_vs,
                                      dual_segment_edge);
      polyscope::registerCurveNetwork("MM dual intersection",
                                      dual_intersect_segement_vs,
                                      dual_intersect_segement_edge);
    }
    medial_mesh->addFaceScalarQuantity("is_dual_seg_intersect",
                                       is_dual_seg_intersect);
    medial_mesh->addFaceScalarQuantity("is_deleted", mat_faces_deleted)
        ->setEnabled(true);
    medial_mesh->addVertexScalarQuantity("radius", mat_radius,
                                         polyscope::DataType::MAGNITUDE);
    // enable back face culling (we dont need to worry about face orientation)
    medial_mesh->setBackFacePolicy(polyscope::BackFacePolicy::Identical);
    medial_mesh->addFaceScalarQuantity("importance", mat_faces_importance)
        ->setEnabled(true);
  }

  ///////////////////////////////////
  // others
  if (polyscope::hasSurfaceMesh("Input mesh"))
    polyscope::getSurfaceMesh("Input mesh")->setEnabled(false);
  if (polyscope::hasPointCloud("sample points"))
    polyscope::getPointCloud("sample points")->setEnabled(false);
  if (polyscope::hasSurfaceMesh("RPD"))
    polyscope::getSurfaceMesh("RPD")->setEnabled(false);
}

void GuiWindow::convert_to_show_edges(
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

}  // namespace matfp