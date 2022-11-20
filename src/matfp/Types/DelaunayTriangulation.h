// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include "matfp/Logger.h"
#include "matfp/Types/CommonTypes.h"

namespace matfp {

// Delaunay Vertex info
class DVI {
 public:
  int tag = -1;  // should be mapping seed_points
  Vector3 pos;
  Vector3 normal = Vector3::Zero();
  bool is_on_feature = false;
};

// Delaunay Tetrahedron info
class DTI {
 public:
  // index is matching shape3D.is_outside
  Vector3 center;
  // tag here is mapping all_medial_spheres, not valid_medial_spheres
  int tag = -1;
  std::string identifier;  // "x-x-x-x", each x is the tag of vertice
                           // (DVI::tag), sorting from low to high
  double sq_radius;
  bool is_outside = true;
  double winding_num;
  // bool is_added_by_dt_for_apporaching_se = false;

  // only inside cell contains a inside_id
  // this is used for calculating initial voronoi diagram
  int inside_id;
};

typedef CGAL::Triangulation_vertex_base_with_info_3<DVI, matfp::Kf> Vb_dt;
typedef CGAL::Triangulation_cell_base_with_info_3<DTI, matfp::Kf> Cb_dt;
typedef CGAL::Triangulation_data_structure_3<Vb_dt, Cb_dt> Tds_dt;
typedef CGAL::Delaunay_triangulation_3<Kf, Tds_dt> Dt;

typedef Dt::Vertex_handle Vertex_handle_dt;
typedef Dt::Cell_handle Cell_handle_dt;

typedef Dt::Vertex_iterator Vertex_iterator_dt;
typedef Dt::Edge_iterator Edge_iterator_dt;
typedef Dt::Facet_iterator Facet_iterator_dt;
typedef Dt::Cell_iterator Cell_iterator_dt;

typedef Dt::Finite_cells_iterator Finite_cells_iterator_dt;
typedef Dt::Finite_facets_iterator Finite_facets_iterator_dt;
typedef Dt::Finite_edges_iterator Finite_edges_iterator_dt;
typedef Dt::Finite_vertices_iterator Finite_vertices_iterator_dt;

typedef Dt::Facet_circulator Facet_circulator_dt;
typedef Dt::Cell_circulator Cell_circulator_dt;
typedef Dt::Point Point_dt;

///////////////
class DelaunayTriangulation : public Dt {
 public:
  // DelaunayTriangulation();
  /**
   * \brief DelaunayTriangulation destructor
   */
  ~DelaunayTriangulation(){};

 public:
  ///////////////
  inline void clean() { this->clear(); }

  inline void set_largest_v_tag(size_t _tag) { largest_v_tag = _tag; }
  inline void set_largest_cell_tag(size_t _cell_tag) {
    largest_cell_tag = _cell_tag;
  }

  inline size_t get_largest_v_tag() const { return largest_v_tag; }

  inline size_t get_largest_cell_tag() { return largest_cell_tag; }

  inline void print_info() {
    logger().info("------ Delaunay Triangulation Info ------");
    logger().info("#v: {}", number_of_vertices());
    logger().info("#e: {}", number_of_finite_edges());
    logger().info("#f: {}", number_of_finite_facets());
    logger().info("#c: {}", number_of_finite_cells());
  }

 private:
  // this should = seed_points + 8 bbox
  // but not != dt.number_of_vertices(), becuase we have some deleted points
  size_t largest_v_tag;
  size_t largest_cell_tag;
};

}  // namespace matfp