// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>
#include <CGAL/Regular_triangulation_vertex_base_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <geogram/basic/common.h>
#include <geogram/basic/counted.h>
#include <geogram/basic/packed_arrays.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/delaunay/delaunay.h>
#include <matfp/Logger.h>
#include <matfp/Types/CommonTypes.h>

#include <unordered_set>

namespace matfp {

// Regular Vertex info
class RVI {
 public:
  bool is_outside = true;
  // tag here is mapping valid_medial_spheres, not all_medial_spheres
  int tag = -1;  // tag is matching index of valid_medial_spheres, not DTI::tag
                 // (because we delete&add mat around features)
  int all_tag = -1;  // matching all_medial_spheres
  std::string identifier;
  Vector3 pos;
  double radius = 0.;
  double sq_radius = 0.;
  double winding_number;  // deprecating
  bool is_feature = false;
  bool is_corner = false;
  bool is_on_bbox = false;  // used for generate_RT_for_dual_PD()
  std::unordered_set<int>
      neig_feature_tags;  // storing neighboring feature vertices tags

  int inside_id = -1;  // only for internal spheres
};

// Regular Tetrahedron info
class RTI {
 public:
  int id = -1;  // assigned while iterating RT, not unique for each run
  bool is_infinite = false;
  std::array<int, 4> all_vs_tags;    // all_idx
  std::array<int, 4> valid_vs_tags;  // valid_idx

  // adjacent RFI::id
  // updated once RFI is stored in rt.rt_fs_info
  // NOTE: rt cell base will not be updated!!
  //       only updated in rt.rt_ts_info
  std::set<int> adj_facets;

  Vector3 center;  // weighted circumcenter
  double dist2sf_dual_center = -1;
  double sq_radius = -1;
  bool is_dual_point_outside = false;
  bool is_contain_corner = false;

  void print_info() const {
    logger().debug("------------ RT tet {}", id);
    logger().debug("all_vs_tags: {},valid_vs_tags: {}", all_vs_tags,
                   valid_vs_tags);
    logger().debug(
        "center: ({},{},{}), sq_radius: {}, dist2sf_dual_center: "
        "{}, is_dual_point_outside: {}, is_contain_corner: {}",
        center[0], center[0], center[0], sq_radius, dist2sf_dual_center,
        is_dual_point_outside, is_contain_corner);
  }
};

typedef CGAL::Regular_triangulation_vertex_base_3<Kf> Vb0_rt;
typedef CGAL::Triangulation_vertex_base_with_info_3<RVI, Kf, Vb0_rt> Vb_rt;
typedef CGAL::Regular_triangulation_cell_base_3<Kf> Cb0_rt;
typedef CGAL::Triangulation_cell_base_with_info_3<RTI, Kf, Cb0_rt> Cb_rt;
typedef CGAL::Triangulation_data_structure_3<Vb_rt, Cb_rt> Tds_rt;
typedef CGAL::Regular_triangulation_3<Kf, Tds_rt> Rt;
typedef CGAL::Triangulation_3<Kf, Tds_rt> Tr;

typedef Kf::Point_3 Point_rt;
typedef Rt::Vertex_iterator Vertex_iterator_rt;
typedef Rt::Vertex_handle Vertex_handle_rt;
typedef Rt::Cell_iterator Cell_iterator_rt;
typedef Rt::Cell_handle Cell_handle_rt;
typedef Rt::Cell_circulator Cell_circulator_rt;
typedef Rt::Facet_circulator Facet_circulator_rt;

typedef Rt::Finite_cells_iterator Finite_cells_iterator_rt;
typedef Rt::Finite_facets_iterator Finite_facets_iterator_rt;
typedef Rt::Finite_edges_iterator Finite_edges_iterator_rt;
typedef Rt::Finite_vertices_iterator Finite_vertices_iterator_rt;
typedef Rt::Tetrahedron Tetrahedron_rt;

// RT face info
// cannot integrate with CGAL
class RFI {
 public:
  int id;            // assigned while iterating RT
  int mat_fid = -1;  // assigned when marked as MAT face
  std::array<int, 3> all_vs_tags = {{-1, -1, -1}};
  std::array<int, 3> valid_vs_tags;
  std::array<int, 2> adj_tets;  // shared by two tets

  // Used for generating MAT
  // see generate_RT_dual_info()
  std::array<Vector3, 2> dual_segment = {{Vector3::Zero(), Vector3::Zero()}};
  std::array<bool, 2> is_dual_seg_vs_outside = {{false, false}};
  // two intersection of dual_segment and surface
  // if not exists (is_dual_seg_intersect = false):
  // store two Vector3::Zero();
  std::array<Vector3, 2> dual_intersections = {
      {Vector3::Zero(), Vector3::Zero()}};
  double dist_dual_intersections = -1;  // distance of intersections [no use]
  bool is_dual_seg_intersect = false;

 public:
  void print_info() const {
    logger().debug("------------ RT face {}", id);
    logger().debug("all_vs_tags: {},valid_vs_tags: {}, adj_tets: {}",
                   all_vs_tags, valid_vs_tags, adj_tets);
    logger().debug(
        "segment: ({},{},{}) -> ({},{},{}), is_dual_seg_vs_outside: {} ",
        dual_segment[0][0], dual_segment[0][1], dual_segment[0][2],
        dual_segment[1][0], dual_segment[1][1], dual_segment[1][2],
        is_dual_seg_vs_outside);
    logger().debug(
        "dual_intersection: ({},{},{}) -> ({},{},{}), dist_dual_intersections: "
        "{}, is_dual_seg_intersect: {}",
        dual_intersections[0][0], dual_intersections[0][1],
        dual_intersections[0][2], dual_intersections[1][0],
        dual_intersections[1][1], dual_intersections[1][2],
        dist_dual_intersections, is_dual_seg_intersect);
  }
};

// RT edge info
// cannot integrate with CGAL
class REI {
 public:
  int id;            // assigned while iterating RT
  int mat_eid = -1;  // assigned when marked as MAT face
  std::array<int, 2> all_vs_tags;
  std::array<int, 2> valid_vs_tags;
  std::vector<int> adj_tets;          // shared by >=3 tets
  std::vector<int> adj_facets;        // shared by 2 faces
  std::vector<Vector3> dual_polygon;  // not triangulized
  std::vector<bool> is_dual_poly_vs_outside;
  // 1. any dual vs is inside
  // 2. dual poly has intersect by
  //    calling aabb_wrapper.sf_polygon_intersection()
  bool is_dual_poly_intersect = false;
  std::vector<int> poly_bbox_isect_triangles;

 public:
  void print_info() const {
    logger().debug("------------ RT edge {}", id);
    logger().debug("all_vs_tags: {},valid_vs_tags: {}, adj_tets: {}",
                   all_vs_tags, valid_vs_tags, adj_tets);
    logger().debug("is_dual_poly_vs_outside: {}, is_dual_poly_intersect: {}",
                   is_dual_poly_vs_outside, is_dual_poly_intersect);
    logger().debug("poly_bbox_isect_triangles: {}", poly_bbox_isect_triangles);
  }
};

///////////////
class RegularTriangulationNN : public Rt, public GEO::Counted {
 public:
  // RegularTriangulationNN();
  /**
   * \brief RegularTriangulationNN destructor
   */
  ~RegularTriangulationNN() {};

  std::vector<RTI> rt_ts_info;  // RT tetras
  std::vector<RFI> rt_fs_info;  // RT faces
  std::vector<REI> rt_es_info;  // RT edges

 public:
  ///////////////
  inline void clean() {
    this->clear();
    tag_to_vh.clear();
    nb_vertices = 0;

    rt_fs_info.clear();
    rt_es_info.clear();
    rt_ts_info.clear();
  }

  inline void print_info() {
    logger().info("------ Regular Triangulation Info ------");
    logger().info("nb_vertices: {}", nb_vertices);
    logger().info("#v: {}", number_of_vertices());
    logger().info("#e: {}", number_of_finite_edges());
    logger().info("#f: {}", number_of_finite_facets());
  }
  inline void print_dual_info() {
    logger().info("------ Regular Triangulation Dual Info ------");
    logger().info("#tets dual vs: {}", rt_ts_info.size());
    logger().info("#f dual segs: {}", rt_fs_info.size());
    logger().info("#e dual polys: {}", rt_es_info.size());
  }

  // ATTENTION:
  // we store valid_medial_spheres.size() here
  // is because rt tag range is [0, valid_medial_spheres.size()]
  // and RPD would use this tag for flagging
  inline void set_nb_vertices(const int nb_v) { nb_vertices = nb_v; }

  inline int get_nb_vertices() const { return nb_vertices; }

  inline Vector3 to_geo_vec(const Point_rt& wp) const {
    return Vector3(wp.x(), wp.y(), wp.z());
  };

  inline double tet_radius(const Tetrahedron_rt& tet) const {
    return (to_geo_vec(tet.vertex(0)) - to_geo_vec(CGAL::circumcenter(tet)))
        .norm();
  };

  inline Vertex_handle_rt get_vh(const GEO::index_t tag) const {
    if (tag_to_vh.empty()) {
      log_and_throw("tag_to_vh cannot be empty");
    } else if (tag_to_vh.find(tag) == tag_to_vh.end()) {
      logger().error("tag {} cannot be found at tag_to_vh", tag);
      log_and_throw("tag not match tag_to_vh value");
    } else if (tag != tag_to_vh.at(tag)->info().tag) {
      logger().error("tag {} has tag_to_vh {}", tag,
                     tag_to_vh.at(tag)->info().tag);
      log_and_throw("tag not match tag_to_vh value");
    } else
      return tag_to_vh.at(tag);
  };

  inline void set_tag_to_vh(int tag, Vertex_handle_rt& vh) {
    tag_to_vh.insert(std::pair<int, Vertex_handle_rt>(tag, vh));
  }

  inline std::vector<double> get_double_vector(const Weighted_point& wp) const {
    std::vector<double> p;
    p.push_back(CGAL::to_double(wp.x()));
    p.push_back(CGAL::to_double(wp.y()));
    p.push_back(CGAL::to_double(wp.z()));
    // this is necessary, since we're calling RPD Vertex::intersect_geom()
    p.push_back(CGAL::to_double(wp.weight()));  // weight,
    return p;
  }

  inline std::vector<double> get_double_vector(const GEO::index_t tag) const {
    const Vertex_handle_rt vh = get_vh(tag);
    return get_double_vector(vh->point());
  };

  /////////////////////////////////////////////////////////
  /////  TODO: try to deprecate following functions
  ///// 		 this type of memory allocation might cause memory leak

  // need to delete[] after calling
  inline double* get_double_data(const Weighted_point& wp) const {
    // new operator creates variable p in the HEAP
    // instead of StACK which willl be destroyed
    // when the function is finished
    // checkout this video
    // https://www.youtube.com/watch?v=RWNM7CzDNyY
    double* p = new double[4];
    p[0] = CGAL::to_double(wp.x());
    p[1] = CGAL::to_double(wp.y());
    p[2] = CGAL::to_double(wp.z());
    // p[3] = 0; // weight
    p[3] = CGAL::to_double(wp.weight());  // weight
    return p;
  };

  // need to delete[] after calling
  inline double* get_double_data(const Vector3& wp) const {
    // new operator creates variable p in the HEAP
    // instead of StACK which willl be destroyed
    // when the function is finished
    // checkout this video
    // https://www.youtube.com/watch?v=RWNM7CzDNyY
    double* p = new double[3];
    p[0] = wp[0];
    p[1] = wp[1];
    p[2] = wp[2];
    return p;
  }

  // need to delete[] after calling
  inline double* get_double_data(const GEO::index_t tag) const {
    Vertex_handle_rt vh = get_vh(tag);
    return get_double_data(vh->point());
  };
  /////////////////////////////////////////////////////////

  inline double get_weight(const Weighted_point& wp) const {
    return CGAL::to_double(wp.weight());
  }
  inline double get_weight(const GEO::index_t tag) const {
    const Vertex_handle_rt& vh = get_vh(tag);
    return CGAL::to_double(vh->point().weight());
  }

 protected:
  // vertex tag -> vertex handle
  std::map<int, Vertex_handle_rt> tag_to_vh;

  // nb_vertices >= number_of_vertices()
  int nb_vertices;
};

/**
 * \brief Smart pointer that refers to a Regular Triangulation object
 */
typedef GEO::SmartPointer<RegularTriangulationNN> RegularTriangulationNN_var;

}  // namespace matfp