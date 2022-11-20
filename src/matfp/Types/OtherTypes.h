// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <set>
#include <unordered_set>
#include <vector>

#include "matfp/Common.h"
#include "matfp/Logger.h"
#include "matfp/Types/CommonTypes.h"

namespace matfp {
/////////////////////////////////////////////////////////////////////////////////////
// TangentConcaveLine
/////////////////////////////////////////////////////////////////////////////////////
// from SpecialEdge with type EdgeType::CCE
// initialized through function load_concave_edges() in MeshFeatureProcessor.h
class TangentConcaveLine {
 public:
  TangentConcaveLine(){};
  ~TangentConcaveLine(){};

  TangentConcaveLine(const int _id, const std::array<Vector3, 2> _end_pts,
                     const std::array<Vector3, 2>& _ns,
                     const std::array<int, 2>& _fpair,
                     const std::array<int, 2>& _vpair);

 public:
  void print_info() const;
  // if concave line is a curve
  // then each line should cover more adjacent faces
  bool is_normal_covered_by_adj_fs(
      const Vector3& n, double esp_degree_given = esp_degree_5) const;
  // given two points of a segment
  void update_tangent_point(const GEO::vec3& v1_p, const GEO::vec3& v2_p);

  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha3, const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha3);

 public:
  bool operator==(TangentConcaveLine const& b) const;

 public:
  int id;
  bool is_tan_point_updated;  // point M will be updated to tangent point X
                              // after each RPD

  Vector3 direction;  // (ref_vs_pos[1] - ref_vs_pos[0]).normalize()
  Vector3 tan_point;  // tangent point, init as center of ref_vs_pos
  // normals is a random vector inside the range of two
  // adjacent normals
  Vector3 normal;

  // distance from sphere to concave line
  // updated by alpha_3
  double energy;
  double energy_over_sq_radius;  // used as break threshold

  // facet id of input GEO::Mesh for RPD
  std::array<int, 2> adj_ref_fs_pair;
  // two normals of adjacent reference plane
  std::array<Vector3, 2> adj_ref_normals;
  // a pair of two vertices of orignal input surface
  // only used for shriking spheres
  std::array<int, 2> ref_vs_seed_ids;
  // start and end point of orignal input surface
  std::array<Vector3, 2> ref_vs_pos;
};

class TangentPlane {
 public:
  TangentPlane();
  TangentPlane(const Vector3& _normal);
  TangentPlane(const Vector3& _normal, const Vector3& _point);
  ~TangentPlane(){};

 public:
  bool clear();
  bool is_same_normal(const Vector3& bnormal,
                      const double eps_degree = esp_degree_10) const;
  static bool is_same_normal(const Vector3& anormal, const Vector3& bnormal,
                             const double eps_degree = esp_degree_10);
  static void union_two_tan_pls_vectors(
      const std::vector<TangentPlane>& tan_pls1,
      const std::vector<TangentPlane>& tan_pls2,
      std::vector<TangentPlane>& tan_pls_merged,
      const double eps_degree = esp_degree_10);
  static void diff_two_tan_pls_vectors(
      const std::vector<TangentPlane>& tan_pls1,
      const std::vector<TangentPlane>& tan_pls2,
      std::vector<TangentPlane>& tan_pls_diff,
      const double eps_degree = esp_degree_10);
  bool is_same_tan_pl(const TangentPlane& tan_pl2) const;
  void push_new_point(Vector3 _p);
  void print_info() const;
  void update_points_with_centroid();
  void update_by_sf_mesh(const GEO::Mesh& sf_mesh,
                         const AABBWrapper& aabb_wrapper);
  static double get_energy_value(const Vector3& theta, const double& radius,
                                 const double alpha1, const double alpha2,
                                 const Vector3& tan_point,
                                 const Vector3& normal);
  double update_energy_value(const Vector3& theta, const double& radius,
                             const double alpha1, const double alpha2);

  void copy_points_from(const TangentPlane& b);
  void sanity_check() const;

  void update_extrude_ratio(const Vector3& center, const double& sq_radius);

 public:
  bool operator==(TangentPlane const& b) const;
  /**
   * Compares two 4D points with respect to the lexicographic order.
   * s1, s2 are the coordinates of the two 4D points.
   * true if s1 is strictly before s2 in the lexicographic order.
   * false otherwise.
   */
  bool operator<(const TangentPlane& t2) const;

 public:
  Vector3 normal;               // normal of tangent plane
  std::vector<Vector3> points;  // points on plane
  double fid;                   // corresponding fid from sf_mesh, -1 as default
  double energy;                // sphere to minimize this energy function,
                                // DBL_MAX default
  double energy_over_sq_radius;  // energy / sq_radius, used as break threshold

  bool is_deleted;
  bool is_touched;

  /*[deprecating]*/
  double extrude_ratio;         // extrude distance / radius, measure how much a
                                // sphere extrudes from the tangent plane
  double extrude_ratio_signed;  // + extrude, - not touching /*[deprecating]*/
};

class ConnectedComponent {
 public:
  ConnectedComponent() {
    is_deleted = false;
    is_on_cc_line = false;
  };
  ~ConnectedComponent(){};

 public:
  void clear();
  void clear_set_storage();

  void print_info() const;
  bool is_adj_to(const int adj_tag) const;

 public:
  std::set<int> adj_sphere_tags;
  std::set<int>
      rpd_facet_indices;  // rpd facet ids that belongs to this cc_part
  std::set<int> rpd_boundary_vertices;  // each RPD vertex has >= 3 bisectors,
                                        // we call it 'boundary'
  bool is_deleted;

  // checking if this cc_part corss a concave line
  bool is_on_cc_line;
  // input facet ids (rpd reference facet ids) that belongs to this cc_part,
  // used for checking normal variance of this cc_part
  std::set<int> rpd_ref_fs_indices;
};

// for ShrinkSpheres
struct ss_params {
  // matching seed_points
  // for making sure different surface sample is selected
  // mostly for unique-ness, no other real purpose
  int pidx = -1;

  // these params are used for real shrinking algo
  Vector3 p, pnormal, q, qnormal;
  int q_fid = -1;       // matching GEO::Mesh facets
  int q_fid_prev = -1;  // matching GEO::Mesh facets
};

// MAT vertices: Data Transfer Object
class MVertex {
 public:
  MVertex(int t, double v0, double v1, double v2, double sq_r = 0.,
          SphereType _st = SphereType::T_UNK);

  MVertex(int t, double v0, double v1, double sq_r = 0.);

  // ShrinkSpheres
  // Init a sphere with relatively large radius
  // and center is not specified
  MVertex(const int t, const Vector3& pin, const Vector3& pin_normal);

 public:
  /**
   * Compares two 4D points with respect to the lexicographic order.
   * s1, s2 are the coordinates of the two 4D points.
   * true if s1 is strictly before s2 in the lexicographic order.
   * false otherwise.
   */
  // SCALAR_ZERO seems to small, use SCALAR_ZERO_N3 instead
  bool operator<(const MVertex& m2) const;
  // for removing medial spheres that are too close
  // using absolute value is not reliable for all models
  // let's use ratio to measure how deep two spheres intersect
  bool operator==(const MVertex& m2) const;
  // for sanity check
  bool operator!=(const MVertex& m2) const;

 public:
  void print_info() const;
  void print_tan_planes() const;
  void print_cc_parts() const;
  void print_cc_lines() const;

  bool ensure_valid();

  int get_num_cc_parts_exclude_cc_lines() const;
  bool is_a_feature_sphere() const;
  bool is_a_concave_sphere() const;
  // if pass concave lines, then add two adjacnet faces' normals
  // and add all tangent planes normals
  void get_sphere_all_tangent_normals_includes_cc_lines(
      std::vector<Vector3>& normals) const;
  void get_sphere_all_tangent_pairs_includes_cc_lines(
      std::vector<std::array<Vector3, 2>>& tan_pairs) const;
  void add_tan_planes_for_feature_sphere(
      const std::vector<Vector3>& normals,
      const std::vector<int>& ref_fids /* matching GEO::Mesh*/);

  bool is_normal_covered_by_tan_pls(
      const Vector3& n, const double esp_degree_given = esp_degree_10) const;
  bool is_normal_covered_by_cc_lines(
      const Vector3& n, const double esp_degree_given = esp_degree_10) const;
  bool is_normal_covered_then_merge(
      const TangentPlane& tan_pl_new,
      const double esp_degree_given = esp_degree_10);
  bool is_normal_covered_includes_cc_lines(
      const Vector3& n, const double esp_degree_given = esp_degree_10) const;
  bool is_normal_covered_by_cc_lines_adj_fs(
      const Vector3& n, const double esp_degree_given) const;
  void update_tan_planes_points();
  void update_tan_planes_by_sf_mesh(const GEO::Mesh& sf_mesh,
                                    const AABBWrapper& aabb_wrapper);

  // for sphere shrinking
  bool update_sphere_ss_param(bool is_debug = false);
  // Initialize sphere to be tangent to pin point p with normal pnormal
  bool init_ss_sphere(bool is_debug = false);
  // update mat_p.tan_planes by p and q from shrinking sphere
  void update_tan_planes_by_ss(bool is_debug = false);

  void save_old_pos_radius(bool is_clear = true);

  double dist_to_plane_signed(const Vector3& pl_point,
                              const Vector3& pl_normal);
  double get_extrude_ratio_unsigned(const Vector3& pl_point,
                                    const Vector3& pl_normal);
  double get_extrude_ratio_signed(const Vector3& pl_point,
                                  const Vector3& pl_normal);
  void dilate_sphere_radius();
  double max_dist_to_tan_points(const std::vector<TangentPlane>& tan_pls);
  double update_all_energy_values(const double alpha1, const double alpha2,
                                  const double alpha3);

 public:
  // tag should be the same as all_medial_spheres index
  // not index of valid_medial_spheres
  int tag;  // unique but not repeatible, same as all_idx
  int valid_idx;
  SphereType type = SphereType::T_UNK;

  // only mat point has identifier
  std::string identifier = "";  // unique and repeatible

  // Local Feature Size
  // 4 numbers that consists of identifier
  // indices are mapping to lfs_points & dt_vs
  // used by estimating LFS
  //
  // Note: this is also used for drawing dt tet with features
  int tet_type = -1;
  std::vector<int> dt_seed_ids;
  std::vector<Vector3> dt_seed_pos;      // ShrinkSpheres
  std::vector<Vector3> dt_seed_normals;  // ShrinkSpheres
  ss_params ss;                          // ShrinkSpheres

  // store all co-incidental Voronoi vertices (mat vs)
  // sphere consists of >5 co-sphere dt seeds
  // including tag of current sphere (size always > 0)
  std::vector<int> co_incident_vs;
  int co_min_tag;  // minimun tag among co_incident_vs
  // Vector4 coord_4d; // pos + radius, use for grouping

  bool is_contains_feature_sample = false;
  size_t num_of_feature_samples = 0;  // [0,4]

  bool is_locking = false;

  ///////////////////////////////////////
  // Update MAT based on tangent planes
  //
  // including both tet and RPD
  std::vector<ConnectedComponent> cc_parts;
  std::vector<TangentPlane> tan_planes;
  std::vector<int> selected_tan_pls;  // a vector of id in tan_planes selected
                                      // for updating this sphere
  std::vector<TangentConcaveLine> tan_cc_lines;
  std::vector<int> selected_cc_lines;
  bool is_to_be_updated = false;   // mark as going to be updated
  bool is_sphere_updated = false;  // will be assign T/F during sphere update
  Vector3 pAB;  // calculated by function get_inter_seg_and_bisector()

  /////////////////////////////////////////
  // To solve unthin
  bool is_create_for_unthin = false;

  ///////////
  // Sharp Edge
  bool is_on_s_edge = false;
  std::unordered_set<int> se_adj_se;  // indices maps to all_medial_spheres
  std::vector<Vector3> se_adj_faces_normals;
  // for debug
  std::set<std::array<int, 2>> k_nearest_se;  // store k nearest se
  std::set<int>
      newly_added_feature_mat_tags;  // store newly added feature points if not
                                     // been deleted, see
                                     // MedialMeshGenerator::add_new_point_to_se()

  //////////
  // Concave edge
  bool is_on_cc_edge = false;  // no use

  //////////
  // Corner
  bool is_on_corner = false;
  bool is_added_for_corner = false;  // for solving corner connectivity
  // relates the sphere to two other spheres (one can be corner sphere)
  // to form a mat face (slab), mapping to MVertex::tag.
  std::set<std::array<int, 2>> added_for_two_spheres;
  // surface point only
  std::vector<int> corner_neighbors_sorted;
  // define a small area around corner, updated in
  // func load_corners_min_length(), same as SpecialCorner
  double min_corner_len = DBL_MAX;

  //////////
  // In/Out
  bool is_outside = false;
  double winding_num = 1.;  // inside
  bool is_on_bbox = false;  // if true then is_outside must be true

  ////////
  // Internal feature edge
  bool is_on_intf_edge = false;
  std::unordered_set<int> intf_adj_intf;   // indices maps to all_medial_spheres
  std::unordered_set<int> k_nearest_nonf;  // all_idx, for less search operation

  // pre-processing before RPDs

  // would override everything
  // valid mat points = !is_outside && !is_deleted
  // deleted for whatever reason, including is_deleted_for_se
  bool is_deleted = false;
  bool is_deleted_for_se = false;  // for solving sharp edges connectivity

  // -1 - unknown
  // 0 - not deleted
  // 1 - deleted for FP external deletion
  // 3 - deleted for cleaning
  DeletionType is_deleted_int =
      DeletionType::D_UNK;  // deleted for whatever reason, including
                            // is_deleted_for_se

  // 0 - not sampled
  // 1 - sampled for FP external corner
  // 2 - sampled for FP internal addition
  AdditionType is_addition_int = AdditionType::A_NOT;

  // for solving two neighboring sharp edges
  bool is_added_by_dt_for_apporaching_se = false;
  int seed_added_idex;  // mapping shap3D.seed_added_for_two_approaching_ses
  bool is_added_for_approaching_se = false;  // non-feature mat point for now
  bool is_added_on_se_for_approaching_se = false;  // feature mat point

  // when is_added_for_approaching_se = true, set these for debug
  // store corresponding 2 feature mat info, and 1 midpoint info
  std::vector<Vector3>
      corres_points;  // 2 feature mat point + 1 projected midpoint info
  std::vector<std::vector<Vector3>>
      corres_normals;  // size matching corres_points

  Vector2 pos_2d = Vector2::Zero();
  Vector3 pos = Vector3::Zero();
  double sq_radius = 0.;
  // dilated radius, increase the radius a bit
  // to solve the degeneration in RPD
  double sq_radius_dilated = 0.;

  Vector3 normal;

  // for changing circum-spheres to inscribed spheres
  // only for tet that contains features samples, for now
  Vector3 old_pos;
  double old_sq_radius;

  ////////
  // for tracking energy for iteration
  // energy = (std::sqrt(sq_distance to nearest surface point) -
  // std::sqrt(sq_radius))^2
  std::vector<double> step_2_sq_energy;
  void store_sq_energy(const AABBWrapper& aabb_wrapper);

  ///////////
  // for optimization
  // NO USE FOR NOW
  // TODO: give initial value
  // pos after initial optimization
  Vector3 new_pos = Vector3::Zero();
  double new_sq_radius = 0.;
  // pin during optimization
  Vector3 pin_pos = Vector3::Zero();
  Vector3 pin_normal = Vector3::Zero();
  double h_sq_score = DBL_MAX;

  // assign a color
  std::array<double, 3> color_rgb = {{((double)rand() / (RAND_MAX)),
                                      ((double)rand() / (RAND_MAX)),
                                      ((double)rand() / (RAND_MAX))}};
};

/////////////////////////////////////////////////////////////////////////////////////
// Corner Preservation
/////////////////////////////////////////////////////////////////////////////////////
//
// helper class
// each sharp edge corresponding a sheet
// part of the sheet could be
class SHEET_SEG {
 public:
  void print_info() const;

 public:
  // normally a sheet has 2 tangent planes
  std::vector<TangentPlane> tan_pls;
  // mapping to ThreeDimensionalShape::tan_cc_lines::id
  std::set<int> tan_cc_line_ids;
};

// Corner Preservation
// Corner MAT structure:
// triangle (corner, node, node_son)
//
// We build a tree structure:
// - Node here indicates the feature sphere (type 1) or seam sphere (type 2).
// - A sheet grows either from sharp edge (node of type 1) or from seam on MAT
// (node of type 2).
// - A sheet consists of multiple sheet segments.
// - A sheet segment has 2 tangent elements (plane or concave point).
class SHEET_NODE {
 public:
  SHEET_NODE(int _id, int _t, int _tag, int _c_tag)
      : id(_id), type(_t), related_mat_tag(_tag), corner_tag(_c_tag) {}
  ~SHEET_NODE() {}
  void print_info() const;

 public:
  int id = -1;
  // 1: sheet extends from sharp edge (leaf node)
  // 2: sheet extends from seam
  int type = -1;
  // if type == 1, sharp feature sphere
  // if type == 2, newly created seam sphere (with is_added_for_corner = true)
  int related_mat_tag = -1;
  // each sheet node must relates to a corner sphere
  int corner_tag = -1;
  // if type == 1
  int sp_edge_id = -1;

  // each sheet segment is defined by 2 tangent elements
  // (tangent planes or tangent concave lines)
  std::vector<SHEET_SEG> sheet_segs;

  // no need to store children
  // we store everything in MVertex::added_for_two_spheres
  //   std::vector<const SHEET_NODE*> children;
};

void export_medial_spheres(const std::string& maname,
                           const std::vector<MVertex>& all_medial_spheres);

// saving step_2_sq_energy
void export_sphere_energy_spline(
    const std::string& maname, const std::vector<MVertex>& all_medial_spheres);

}  // namespace matfp
