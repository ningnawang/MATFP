/*
 *  Copyright (c) 2012-2014, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     Bruno.Levy@inria.fr
 *     http://www.loria.fr/~levy
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 */

#pragma once

#include <geogram/basic/argused.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/common.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/process.h>
#include <geogram/mesh/index.h>
#include <geogram/numerics/predicates.h>

#include <algorithm>
#include <deque>
#include <iostream>

#include "matfp/Logger.h"
#include "matfp/Types/RegularTriangulation.h"
#include "matfp/geogram/RPD_callback.h"
#include "matfp/geogram/generic_RPD_cell.h"
#include "matfp/geogram/generic_RPD_polygon.h"
#include "matfp/geogram/generic_RPD_utils.h"
#include "matfp/geogram/generic_RPD_vertex.h"

namespace matfp {
using namespace GEOGen;

/**
 * \brief Symbolic representation of a RestrictedVoronoiDiagram vertex.
 */
typedef matfp::SymbolicVertex SymbolicVertex;

class GenRestrictedPowerDiagram {
  /** \brief This class type */
  typedef GenRestrictedPowerDiagram thisclass;

 public:
  /**
   * \brief Gets the dimension
   */
  static coord_index_t dimension() { return 3; }

  /**
   * \brief Used to allocate the generated points.
   */
  typedef GEOGen::PointAllocator PointAllocator;

  /**
   * \brief Internal representation of vertices.
   */
  typedef matfp::Vertex Vertex;

  /**
   * \brief Internal representation of polygons.
   */
  typedef matfp::PolygonCGAL Polygon;

  /**
   * \brief Internal representation of volumetric cells.
   */
  typedef matfp::ConvexCellCGAL Polyhedron;

  /********************************************************************/

  /**
   * \brief Constructs a new GenRestrictedPowerDiagram.
   * \param[in] rt the regular triangulation
   * \param[in] mesh the input mesh
   */
  GenRestrictedPowerDiagram(RegularTriangulationNN* rt, GEO::Mesh* mesh)
      : mesh_(mesh),
        rt_(rt),
        intersections_(3),
        symbolic_(false),
        check_SR_(true),
        exact_(false) {
    dimension_ = 3;  // though we have weight, dimension should still be 3
    facets_begin_ = UNSPECIFIED_RANGE;
    facets_end_ = UNSPECIFIED_RANGE;
    tets_begin_ = UNSPECIFIED_RANGE;
    tets_end_ = UNSPECIFIED_RANGE;
    connected_components_priority_ = false;
    // facet_seed_marking_ = nullptr;
    // connected_component_changed_ = false;
    // current_connected_component_ = 0;
    // cur_stamp_ = -1;
    current_facet_ = GEO::max_index_t();
    current_seed_ = GEO::max_index_t();
    current_seed_handle_ = nullptr;
    current_polygon_ = nullptr;
    // current_tet_ = GEO::max_index_t();
    // current_polyhedron_ = nullptr;
  }

  /**
   * \brief Sets traveral priority.
   * \details If connected_components_priority is set,
   * then the connected components of the
   * restricted Voronoi cells will be traversed
   * one by one.
   */
  void set_connected_components_priority(bool x) {
    connected_components_priority_ = x;
  }

  /**
   * \brief Tests whether connected components priority is
   *  set.
   * \details If connected_components_priority is set,
   *  then the connected components of the
   *  restricted Voronoi cells will be traversed
   *  one by one.
   * \retval true if connected components priority is used.
   * \retval false otherwise.
   */
  bool connected_components_priority() const {
    return connected_components_priority_;
  }

  /**
   * \brief Gets the index of the mesh facet currently processed.
   * \details Can be used in surfacic traversals (and not volumetric
   *  traversals).
   */
  GEO::index_t current_facet() const { return current_facet_; }
  /**
   * \brief Gets the index of the Delaunay vertex currently processed.
   * \details Can be used in both surfacic traversals and volumetric
   *  traversals.
   */
  GEO::index_t current_seed() const { return current_seed_; }

  Vertex_handle_rt current_seed_handle() const { return current_seed_handle_; }

  /**
   * \brief Gets the current polygon.
   * \details The current polygon corresponds to the
   *  intersection between the current facet
   *  and the Voronoi cell of the current seed. Can be used
   *  in surfacic traversals (and not volumetric traversals).
   */
  const Polygon& current_polygon() const { return *current_polygon_; }

  /********************************************************************/

 public:
  /**
   * @}
   * \name Public interface for computation/iteration
   * @{
   */

  /**
   * \brief Iterates on the facets of this RPD.
   * \param[in] action the user action object
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, index_t f, const Polygon& P) const
   *  where v denotes the index of the current Voronoi cell
   *  (or Delaunay vertex), f the index of the current facet
   *  and P the computed intersection between facet f
   *  and the Voronoi cell of v.
   */
  template <class ACTION>
  inline void for_each_polygon(const ACTION& action) {
    // PolygonAction<ACTION> polyaction(action);
    // std::cout << "+++ calling compute_surfacic_with_seeds_priority ... \n";
    this->template compute_surfacic_with_seeds_priority<PolygonAction<ACTION>>(
        PolygonAction<ACTION>(action));
    // std::cout << "finish for_each_polygon" << std::endl;
  }

  /**
   * \brief Iterates on the facets of this RPD, triangulated on the fly.
   * \param[in] action the user action object
   * \tparam TRIACTION needs to implement:
   *  operator()(index_t c, const TopoPolyVertex& v1, v2, v3) const
   *  where c denotes the index of the current Voronoi cell
   *  (or Delaunay vertex).
   */
  template <class TRIACTION>
  inline void for_each_triangle(const TRIACTION& action) {
    this->template compute_surfacic_with_seeds_priority<
        TriangleAction<TRIACTION>>(TriangleAction<TRIACTION>(action));
    // std::cout << "finish for_each_triangle" << std::endl;
  }

  /**
   * \brief Iterates on the polyhedra of this RVD.
   * \param[in] action the user action object
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, index_t t, const Polyhedron& C) const
   *  where v denotes the index of the current Voronoi cell
   *  (or Delaunay vertex), t the index of the current tetrahedron
   *  and C the computed intersection between tetrahedron t
   *  and the Voronoi cell of v.
   */
  template <class ACTION>
  inline void for_each_polyhedron(const ACTION& action) {
    this->template compute_volumetric_with_seeds_priority<
        PolyhedronAction<ACTION>>(PolyhedronAction<ACTION>(action));
  }

  /**
   * \brief Iterates on the polyhedra of this RVD decomposed
   *  on the fly into tetrahedra.
   * \details The generated tetrahedra may be geometrically incorrect,
   *  but they are algebraically correct. In other word, their signed
   *  volumes sum as the volume of the restricted Voronoi cell.
   * \param[in] action the user action object
   * \param[in] visit_inner_tets if set, all the tetrahedron-cell
   *  intersections are visited, else only tetrahedra on the border
   *  of the restricted Voronoi cell are visited. Since all the visited
   *  triangles are connected to the current Voronoi seed by a
   *  tetrahedron, the computed volume is the same in both cases.
   * \param[in] coherent_triangles if set, this ensures that the
   *  polygonal facets of the cells are always triangulated in a
   *  coherent manner when seen from two different cells.
   *  For instance, it is required if a tetrahedral mesh is
   *  reconstructed.
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, signed_index_t v_adj,
   *    index_t t, index_t t_adj,
   *    const Vertex& v1, const Vertex& v2, const Vertex& v3
   *  )
   *  where the parameters are as follows:
   *    - v is the index of the current Voronoi cell
   *    (or Delaunay vertex)
   *    - v_adj is the index of the Voronoi cell adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    adjacent to v or -1 if current face is a tetrahedron facet
   *    - t is the index of the current tetrahedron
   *    - t_adj is the index of the tetrahedron adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    - v1,v2 and v3 are the three vertices of the facet on the
   *    border of the restricted Voronoi cell.
   */
  template <class ACTION>
  inline void for_each_volumetric_integration_simplex(
      const ACTION& action, bool visit_inner_tets = false,
      bool coherent_triangles = false) {
    this->template compute_volumetric_with_seeds_priority<
        VolumetricIntegrationSimplexAction<ACTION>>(
        VolumetricIntegrationSimplexAction<ACTION>(action, visit_inner_tets,
                                                   coherent_triangles));
  }

  /**
   * \brief Iterates on the polyhedra of this RPD decomposed
   *  on the fly into tetrahedra.
   * \details The tetrahedra are generated by connecting one of
   *  the vertices of the cell to the other ones.
   * \param[in] action the user action object
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, signed_index_t v_adj,
   *    index_t t, index_t t_adj,
   *    const Vertex& v0, const Vertex& v1,
   *    const Vertex& v2, const Vertex& v3
   *  )
   *  where the parameters are as follows:
   *    - v is the index of the current Voronoi cell
   *    (or Delaunay vertex)
   *    - v_adj is the index of the Voronoi cell adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    adjacent to v or -1 if current face is a tetrahedron facet
   *    - t is the index of the current tetrahedron
   *    - t_adj is the index of the tetrahedron adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    - v0,v1,v2 and v3 are the four vertices of tetrahedron.
   */
  template <class ACTION>
  inline void for_each_tetrahedron(const ACTION& action) {
    this->template compute_volumetric_with_seeds_priority<
        TetrahedronAction<ACTION>>(TetrahedronAction<ACTION>(action));
  }

  template <class ACTION>
  inline void for_given_tet_and_seed(
      const GEO::index_t& tet_idx, const Vertex_handle_rt& current_seed_handle,
      const Vertex_handle_rt& neigh_seed_handle, const ACTION& action,
      const bool integration_simplices) {
    if (integration_simplices) {
      this->template compute_volumetric_given_tet_and_seed<
          VolumetricIntegrationSimplexAction<ACTION>>(
          tet_idx, current_seed_handle, neigh_seed_handle,
          VolumetricIntegrationSimplexAction<ACTION>(action, false, true));
    } else {
      this->template compute_volumetric_given_tet_and_seed<
          TetrahedronAction<ACTION>>(tet_idx, current_seed_handle,
                                     neigh_seed_handle,
                                     TetrahedronAction<ACTION>(action));
    }
  }

  /********************************************************************/

 public:
  /**
   * \brief Low-level API of Restricted Voronoi Diagram traversal
   *  with seeds priority in surfacic mode.
   * \details Client code may use for_each_facet(),for_each_triangle() or
   *  for_each_primal_triangle() instead.
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, index_t f, const Polygon& P) const
   *  where v denotes the index of the current Voronoi cell
   *  (or Delaunay vertex), f the index of the current facet
   *  and P the computed intersection between the Voronoi cell of
   *  v and facet f.
   */
  template <class POLYGONACTION>
  inline void compute_surfacic_with_seeds_priority(
      const POLYGONACTION& polygon_action) {
    if (facets_begin_ == UNSPECIFIED_RANGE &&
        facets_end_ == UNSPECIFIED_RANGE) {
      facets_begin_ = 0;
      facets_end_ = mesh_->facets.nb();
    }

    // std::cout << "rt_->get_nb_vertices(): " << rt_->get_nb_vertices() <<
    // std::endl;

    typename GenRestrictedPowerDiagram::Polygon Facet;
    current_polygon_ = nullptr;
    // seed -> facet
    GEO::vector<index_t> seed_stamp(rt_->get_nb_vertices(), index_t(-1));
    GEO::vector<bool> facet_is_marked(facets_end_ - facets_begin_, false);

    FacetSeedHandleStack adjacent_facets;
    SeedHandleStack adjacent_seeds;
    Polygon F;
    GEO::Attribute<double> vertex_weight;
    vertex_weight.bind_if_is_defined(mesh_->vertices.attributes(), "weight");

    // The algorithm propagates along both the facet-graph of
    // the surface and the 1-skeleton of the Delaunay triangulation,
    // and computes all the relevant intersections between
    // each Voronoi cell and facet.
    for (index_t f = facets_begin_; f < facets_end_; f++) {
      // for(index_t f = 33; f < 34; f++) {
      // GEO::index_t f = 1;
      // logger().debug("going to process facet {}", f);

      if (!facet_is_marked[f - facets_begin_]) {
        // Propagate along the facet-graph.
        facet_is_marked[f - facets_begin_] = true;
        adjacent_facets.push(FacetSeedHandle(f, find_seed_near_facet(f)));
        while (!adjacent_facets.empty()) {
          current_facet_ = adjacent_facets.top().f;
          current_seed_handle_ = adjacent_facets.top().seed;
          current_seed_ = current_seed_handle_->info().tag;
          adjacent_facets.pop();

          // logger().debug("symbolic_: {}", symbolic_);
          // logger().debug("--------------processing current facet {} seed {}",
          //     current_facet_, current_seed_
          // );

          // Copy the current facet from the Mesh into
          // RestrictedPowerDiagram's Polygon data structure
          // (gathers all the necessary information)
          Facet.initialize_from_mesh_facet(mesh_, current_facet_, symbolic_,
                                           vertex_weight);

          // Propagate along the Delaunay 1-skeleton
          // This will traverse all the seeds such that their
          // Voronoi cell has a non-empty intersection with
          // the current facet.
          // std::cout << "seed_stamp size: " << seed_stamp.size() << std::endl;
          seed_stamp[current_seed_] = current_facet_;
          // std::cout << "adjacent_seeds size: " << adjacent_seeds.size() <<
          // std::endl;
          adjacent_seeds.push(current_seed_handle_);

          while (!adjacent_seeds.empty()) {
            // Attention: Once we update current_seed_handle_,
            // we need to update current_seed_ as well.
            current_seed_handle_ = adjacent_seeds.top();
            current_seed_ = current_seed_handle_->info().tag;
            adjacent_seeds.pop();

            // if (f == 33) {
            //     logger().debug("--------------processing current facet {}
            //     seed {}",
            //         current_facet_, current_seed_
            //     );
            // }

            // std::cout << "------------------processing current facet " <<
            // current_facet_
            //     << " and seed " << current_seed_ << std::endl;

            current_polygon_ =
                intersect_cell_facet(current_seed_handle_, Facet);

            // std::cout << "after calling intersect_cell_facet" << std::endl;
            // if (current_polygon_ != nullptr)
            //     std::cout << "current_polygon_ is not null, has nb_vertices:
            //     " <<  current_polygon_->nb_vertices() << std::endl;
            // else
            //     std::cout << "current_polygon_ is not null" << std::endl;

            // logger().debug("DONE current_seed_: {}, current_facet_: {},
            // current_polygon: {}",
            //     current_seed_, current_facet_,
            //     current_polygon().nb_vertices());

            // if (current_seed_ == 457) {
            //     logger().debug("--------------processing current facet {}
            //     seed {}",
            //         current_facet_, current_seed_
            //     );
            //     logger().debug("current_polygon_ nb_vertices {}",
            //     current_polygon_->nb_vertices()); logger().debug("DONE
            //     current_seed_: {}, current_facet_: {}, current_polygon: {}",
            //         current_seed_, current_facet_,
            //         current_polygon().nb_vertices());
            // }

            polygon_action(current_seed_, current_facet_, current_polygon());

            // std::cout << "after calling action, propogate to polygon with "
            // << current_polygon().nb_vertices() << " vertices" << std::endl;

            // Propagate to adjacent facets and adjacent seeds
            for (index_t v = 0; v < current_polygon().nb_vertices(); v++) {
              // std::cout << "propogating to polygon v: " << v << std::endl;
              const Vertex& ve = current_polygon().vertex(v);
              GEO::signed_index_t neigh_f = ve.adjacent_facet();
              // std::cout << "neigh_f: " << neigh_f << std::endl;

              if (neigh_f >= GEO::signed_index_t(facets_begin_) &&
                  neigh_f < GEO::signed_index_t(facets_end_) &&
                  neigh_f != GEO::signed_index_t(current_facet_)) {
                if (!facet_is_marked[index_t(neigh_f) - facets_begin_]) {
                  facet_is_marked[index_t(neigh_f) - facets_begin_] = true;
                  adjacent_facets.push(FacetSeedHandle(GEO::index_t(neigh_f),
                                                       current_seed_handle_));
                  // std::cout << "pushed to adjacent_facets" << std::endl;
                }
              }
              GEO::signed_index_t neigh_s = ve.adjacent_seed();

              // std::cout << "neigh_s " << neigh_s << std::endl;

              if (neigh_s != -1) {
                // std::cout << "seed_stamp[neigh_s]: " << seed_stamp[neigh_s]
                // << " current_facet_: " << current_facet_ << std::endl;
                if (seed_stamp[neigh_s] != current_facet_) {
                  seed_stamp[neigh_s] = current_facet_;
                  adjacent_seeds.push(
                      find_adjacent_seed(GEO::index_t(neigh_s)));
                }
              }
              // std::cout << "propogation done" << std::endl;

            }  // propogate to adjacent facets end
               // std::cout << "finish propogation loop" << std::endl;
          }    // while adjacent_seeds end
               // std::cout << "finish adjacent_seeds loop" << std::endl;
        }      // while adjacent_facets end
               // std::cout << "finish adjacent_facets loop" << std::endl;
      }        // if facet_is_marked
               // std::cout << "finish facet_is_marked" << std::endl;
    }          // for facet loop end

    // std::cout << "finish facet loop" << std::endl;

    current_polygon_ = nullptr;

    // std::cout << "finish compute_surfacic_with_seeds_priority" << std::endl;
  }

  /**
   * \brief Low-level API of Restricted Voronoi Diagram traversal
   *  with seeds priority in volumetric mode.
   * \details Client code may use for_each_polyhedron() or
   *  for_each_volumetric_integration_simplex() instead.
   * \tparam ACTION needs to implement:
   *  operator()(index_t v, index_t t, const Polyhedron& C) const
   *  where v denotes the index of the current Voronoi cell
   *  (or Delaunay vertex), t the index of the current tetrahedron
   *  and C the computed intersection between the Voronoi cell of
   *  v and tetrahedron t
   */
  template <class ACTION>
  inline void compute_volumetric_with_seeds_priority(const ACTION& action) {
    if (tets_begin_ == UNSPECIFIED_RANGE && tets_end_ == UNSPECIFIED_RANGE) {
      tets_begin_ = 0;
      tets_end_ = mesh_->cells.nb();
    }
    geo_assert(tets_begin_ != UNSPECIFIED_RANGE);
    geo_assert(tets_end_ != UNSPECIFIED_RANGE);

    GEO::vector<GEO::index_t> seed_stamp(rt_->get_nb_vertices(),
                                         GEO::index_t(-1));
    GEO::vector<bool> tet_is_marked(tets_end_ - tets_begin_, false);
    // init_get_neighbors();

    TetSeedHandleStack adjacent_tets;
    SeedHandleStack adjacent_seeds;
    Polyhedron C(dimension());
    GEO::Attribute<double> vertex_weight;
    vertex_weight.bind_if_is_defined(mesh_->vertices.attributes(), "weight");

    current_polyhedron_ = &C;
    // The algorithm propagates along both the facet-graph of
    // the surface and the 1-skeleton of the Delaunay triangulation,
    // and computes all the relevant intersections between
    // each Voronoi cell and facet.
    for (GEO::index_t t = tets_begin_; t < tets_end_; ++t) {
      // for(GEO::index_t t = tets_begin_; t < 1; ++t) {
      if (!tet_is_marked[t - tets_begin_]) {
        // Propagate along the tet-graph.
        tet_is_marked[t - tets_begin_] = true;
        adjacent_tets.push(TetSeedHandle(t, find_seed_near_tet(t)));
        while (!adjacent_tets.empty()) {
          current_tet_ = adjacent_tets.top().f;
          current_seed_handle_ = adjacent_tets.top().seed;
          current_seed_ = current_seed_handle_->info().tag;
          adjacent_tets.pop();

          // Note: current cell could be looked up here,
          // (from current_tet_) if we chose to keep it
          // and copy it right before clipping (I am
          // not sure that it is worth it, lookup time
          // will be probably fast enough)

          // Propagate along the Delaunay 1-skeleton
          // This will traverse all the seeds such that their
          // Voronoi cell has a non-empty intersection with
          // the current facet.
          seed_stamp[current_seed_] = current_tet_;
          adjacent_seeds.push(current_seed_handle_);

          while (!adjacent_seeds.empty()) {
            current_seed_handle_ = adjacent_seeds.top();
            current_seed_ = current_seed_handle_->info().tag;
            adjacent_seeds.pop();

            // logger().debug("--------------processing current tet {} seed {}",
            //     current_tet_, current_seed_);

            C.initialize_from_mesh_tetrahedron(mesh_, current_tet_, symbolic_,
                                               vertex_weight);

            intersect_cell_cell(current_seed_handle_, C);

            action(current_seed_, current_tet_, current_polyhedron());

            // Propagate to adjacent tets and adjacent seeds
            // Iterate on the vertices of the cell (remember:
            // the cell is represented in dual form)
            for (GEO::index_t v = 0; v < current_polyhedron().max_v(); ++v) {
              //  Skip clipping planes that are no longer
              // connected to a cell facet.
              if (current_polyhedron().vertex_triangle(v) == -1) {
                continue;
              }

              GEO::signed_index_t id = current_polyhedron().vertex_id(v);
              if (id > 0) {
                // Propagate to adjacent seed
                GEO::index_t neigh_s = GEO::index_t(id - 1);
                if (seed_stamp[neigh_s] != current_tet_) {
                  seed_stamp[neigh_s] = current_tet_;
                  adjacent_seeds.push(find_adjacent_seed(neigh_s));
                }
              } else if (id < 0) {
                // id==0 corresponds to facet on boundary
                //  (skipped)
                // id<0 corresponds to adjacent tet index

                // Propagate to adjacent tet
                GEO::signed_index_t neigh_t = -id - 1;
                if (neigh_t >= GEO::signed_index_t(tets_begin_) &&
                    neigh_t < GEO::signed_index_t(tets_end_) &&
                    neigh_t != GEO::signed_index_t(current_tet_)) {
                  if (!tet_is_marked[GEO::index_t(neigh_t) - tets_begin_]) {
                    tet_is_marked[GEO::index_t(neigh_t) - tets_begin_] = true;
                    adjacent_tets.push(TetSeedHandle(GEO::index_t(neigh_t),
                                                     current_seed_handle_));
                  }
                }
              }  // if done
            }    // propogate to adjacent tet done
          }      // while adjacent_seeds end
        }        // while adjacent_tets end
      }          // if tet_is_marked
    }            // for tet loop end
    current_polyhedron_ = nullptr;
  }

  template <class ACTION>
  inline void compute_volumetric_given_tet_and_seed(
      const GEO::index_t& tet_idx, const Vertex_handle_rt& current_seed_handle,
      const Vertex_handle_rt& neigh_seed_handle, const ACTION& action) {
    TetSeedHandleStack adjacent_tets;
    SeedHandleStack adjacent_seeds;
    Polyhedron C(dimension());
    GEO::Attribute<double> vertex_weight;
    vertex_weight.bind_if_is_defined(mesh_->vertices.attributes(), "weight");

    current_polyhedron_ = &C;
    // The algorithm propagates along both the facet-graph of
    // the surface and the 1-skeleton of the Delaunay triangulation,
    // and computes all the relevant intersections between
    // each Voronoi cell and facet.
    GEO::index_t t = tet_idx;

    adjacent_tets.push(TetSeedHandle(t, current_seed_handle));
    while (!adjacent_tets.empty()) {
      current_tet_ = adjacent_tets.top().f;
      current_seed_handle_ = adjacent_tets.top().seed;
      current_seed_ = current_seed_handle_->info().tag;
      adjacent_tets.pop();

      // Note: current cell could be looked up here,
      // (from current_tet_) if we chose to keep it
      // and copy it right before clipping (I am
      // not sure that it is worth it, lookup time
      // will be probably fast enough)

      // Propagate along the Delaunay 1-skeleton
      // This will traverse all the seeds such that their
      // Voronoi cell has a non-empty intersection with
      // the current facet.
      // seed_stamp[current_seed_] = current_tet_;
      adjacent_seeds.push(current_seed_handle_);

      while (!adjacent_seeds.empty()) {
        current_seed_handle_ = adjacent_seeds.top();
        current_seed_ = current_seed_handle_->info().tag;
        adjacent_seeds.pop();

        logger().debug("--------------processing current tet {} seed {}",
                       current_tet_, current_seed_);

        C.initialize_from_mesh_tetrahedron(mesh_, current_tet_, symbolic_,
                                           vertex_weight);

        if (neigh_seed_handle == nullptr) {
          intersect_cell_cell(current_seed_handle_, C);
        } else {
          Vertex_handle_rt neigh_seed_handle_ = neigh_seed_handle;
          clip_by_cell_given_neigh(current_seed_handle_, neigh_seed_handle_, C);
        }
        action(current_seed_, current_tet_, current_polyhedron());

      }  // while adjacent_seeds end
    }    // while adjacent_tets end

    current_polyhedron_ = nullptr;
  }

  /********************************************************************/

  /**
   * @}
   * \name Clipping for surfacic mode
   * @{
   */

 public:
  /**
   * \brief Computes the intersection between the Voronoi cell
   * of a seed and a facet.
   * \param[in] seed the index of the seed
   * \param[in] F the facet represented as a Polygon
   * \details The result is provided in current_polygon_
   */
  Polygon* intersect_cell_facet(Vertex_handle_rt& seed_handle, Polygon& F) {
    intersections_.clear();

    // Initialize ping-pong pointers for Sutherland-Hodgman
    // re-entrant clipping and copy current facet into 'ping' buffer.
    Polygon* ping = &F;
    Polygon* pong = &P2;

    // Clip current facet by current Voronoi cell (associated with seed)
    // if(delaunay_nn_ != nullptr) {
    //     clip_by_cell_SR(seed, ping, pong);   // "Security Radius" mode.
    // } else {
    clip_by_cell(seed_handle, ping, pong);  // Standard mode.
    // }

    return ping;  // Yes, 'ping', and not 'pong'
                  // see comments in clip_by_cell()
  }

 protected:
  /**
   * \brief Swaps two pointers between two polygons.
   * \details Used by re-entrant Sutherlang-Hogdman clipping.
   */
  void swap_polygons(Polygon*& ping, Polygon*& pong) {
    if (ping != &P1 && ping != &P2) {
      // First clipping operation, ping points to F
      // (current facet copied)
      ping = &P2;
      pong = &P1;
    } else {
      // logger().debug("std::swapping");
      std::swap(ping, pong);
    }
  }

  /**
   * \brief Computes the intersection between a Voronoi cell
   *  and a polygon.
   *
   * \details The Voronoi cell is determined by vertex \p i and
   * the input polygon is in \p ping. The result is returned
   *  in \p ping (Note that
   * \p ping and \p pong are references, and that they are swapped
   * after each bisector clipping, this is why the final result
   * is in \p ping (and not in \p pong).
   *
   * \param[in] i index of the vertex that defines the Voronoi cell
   * \param[in,out] ping the input polygon. On exit, contains the result.
   * \param[out] pong a buffer used to implement reentrant clipping.
   *  Its content is modified by the function.
   */
  void clip_by_cell(Vertex_handle_rt& i, Polygon*& ping, Polygon*& pong) {
    get_neighbors(i);
    // logger().debug("seed {} has neighbors_ size {}", i->info().tag,
    // neighbors_.size()); std::cout << "seed " << i->info().tag << " has
    // neighbors size " << neighbors_.size() << std::endl;
    for (index_t jj = 0; jj < neighbors_.size(); jj++) {
      Vertex_handle_rt& j = neighbors_[jj];
      clip_by_plane(*ping, *pong, i, j);
      swap_polygons(ping, pong);
    }
  }

  /**
   * \brief Computes the intersection between a polygon and a half-space.
   *
   * \details The input polygon is in \p ping
   * and the half-space is determined by the positive side
   * of the bisector of segment [\p i,\p j] (the side of \p i).
   * The result is stored into the Polygon \p pong.
   *
   * \param[in] i index of the first extremity of the bisector
   * \param[in] j index of the second extremity of the bisector
   * \param[in] ping the input polygon
   * \param[out] pong \p ping clipped by the bisector
   */

  void clip_by_plane(Polygon& ping, Polygon& pong, Vertex_handle_rt& i,
                     Vertex_handle_rt& j) {
    ping.clip_by_plane<3>(pong, intersections_, mesh_, rt_, i, j, exact_,
                          symbolic_);
  }

  /********************************************************************/

  /**
   * @}
   * \name Clipping for volumetric mode
   * @{
   */

 public:
  /**
   * \brief Computes the intersection between a Voronoi cell
   *  and a cell with radius of security or plain mode.
   * \param[in] seed the index of the seed that defines the Voronoi cell
   * \param[in,out] C the cell to be clipped
   */
  void intersect_cell_cell(Vertex_handle_rt& seed, Polyhedron& C) {
    // Clip current facet by current Voronoi cell (associated with seed)
    // if(delaunay_nn_ != nullptr) {
    // clip_by_cell_SR(seed, C);   // "Security Radius" mode.
    // } else {
    clip_by_cell(seed, C);  // Standard mode.
                            // }
  }

 protected:
  /**
   * \brief Computes the intersection between a Voronoi cell
   *  and a cell in plain mode.
   * \param[in] seed the index of the seed that defines the Voronoi cell
   * \param[in,out] C the cell to be clipped
   */
  void clip_by_cell(Vertex_handle_rt& seed, Polyhedron& C) {
    get_neighbors(seed);

    // Check whether cell is empty (may happen with
    // power diagrams)
    if (neighbors_.size() == 0) {
      C.clear();
    }

    std::vector<int> neighbors;
    for (index_t jj = 0; jj < neighbors_.size(); jj++) {
      Vertex_handle_rt& j = neighbors_[jj];
      neighbors.push_back(j->info().tag);
    }
    // logger().debug("tag {} has {} neighbors {}", seed->info().all_tag,
    // neighbors_.size(), neighbors);

    for (index_t jj = 0; jj < neighbors_.size(); jj++) {
      Vertex_handle_rt& j = neighbors_[jj];
      // logger().debug("processing clip_by_plane i {} and j {}",
      // seed->info().all_tag, j->info().all_tag);
      clip_by_plane(C, seed, j);
    }
    // logger().debug("tag {} has {} neighbors {}", seed->info().all_tag,
    // neighbors_.size(), neighbors);
  }

  void clip_by_cell_given_neigh(Vertex_handle_rt& seed,
                                Vertex_handle_rt& neigh_seed, Polyhedron& C) {
    get_neighbors(seed);

    // Check whether cell is empty (may happen with
    // power diagrams)
    if (neighbors_.size() == 0) {
      C.clear();
    }

    logger().debug("processing clip_by_plane i {} and j {}",
                   seed->info().all_tag, neigh_seed->info().all_tag);
    clip_by_plane(C, seed, neigh_seed);
  }

  /**
   * \brief Computes the intersection between a Voronoi cell
   *  and a half-space determined by a bisector.
   * \param[in,out] C cell to be clipped
   * \param[in] i index of the first extremity of the bisector
   * \param[in] j index of the second extremity of the bisector
   */
  void clip_by_plane(Polyhedron& C, Vertex_handle_rt& i, Vertex_handle_rt& j) {
    C.clip_by_plane<3>(mesh_, rt_, i, j, exact_, symbolic_);
  }

  /********************************************************************/
 public:
  /**
   * @}
   * \name Optimized get neighbors
   * @{
   */

  /**
   * \brief Caches the neighbors of a Delaunay vertex.
   *
   * \details This function is only used when the stored delaunay
   *  triangulation is a traditional one. When the stored delaunay
   *  triangulation is represented by a KdTree, function is not used.
   */
  void get_neighbors(Vertex_handle_rt& v) {
    neighbors_.resize(0);
    rt_->finite_adjacent_vertices(v, std::back_inserter(neighbors_));
    // sort from small to big tag
    std::sort(neighbors_.begin(), neighbors_.end(),
              [](Vertex_handle_rt& a, Vertex_handle_rt& b) {
                return a->info().tag > b->info().tag;
              });
  }

  /********************************************************************/
 protected:
  /**
   * \brief Fetch neighbor seed handle given tag
   *
   * \details This function can only be run under the condition that
   * all neighbors of a given vertex handle are cached
   **/
  Vertex_handle_rt find_adjacent_seed(GEO::index_t neigh_s_tag) {
    if (neighbors_.size() == 0) {
      log_and_throw("neighbors_ not been cached properly!");
    }
    // std::cout << "in find_adjacent_seed, neighbors_ size: " <<
    // neighbors_.size() << std::endl;

    // todo; use the map in RT instead of looping here
    for (Vertex_handle_rt n : neighbors_) {
      if (neigh_s_tag == n->info().tag) return n;
    }

    for (Vertex_handle_rt n : neighbors_) {
      std::cout << "neigh_s_tag: " << neigh_s_tag
                << " and neighbors_ tag: " << n->info().tag << std::endl;
    }
    log_and_throw("neighbor vertex handle not found!");
  }

  /**
   * \brief Finds a seed near a given facet.
   * \param[in] f index of the facet in the mesh
   * \return the index of a Voronoi seed such that there is a
   *   non-empty intersection between the Voronoi cell
   *   of the seed and facet \p f.
   */
  Vertex_handle_rt find_seed_near_facet(index_t f) {
    const double* p = mesh_->vertices.point_ptr(mesh_->facets.vertex(f, 0));
    return find_seed_near_point(p);
  }

  /**
   * \brief Finds a seed near a given tetrahedron.
   * \param[in] t index of the tetrahedron in the mesh
   * \return the index of a Voronoi seed such that there is a
   *   non-empty intersection between the Voronoi cell
   *   of the seed and tetrahedron \p t.
   */
  Vertex_handle_rt find_seed_near_tet(index_t t) {
    index_t v = mesh_->cells.tet_vertex(t, 0);
    const double* p = mesh_->vertices.point_ptr(v);
    return find_seed_near_point(p);
  }

  /**
   * \brief Finds a seed near a given point.
   * \param[in] p pointer to the coordinates of the point
   * \return the index of a Voronoi seed such that its
   *   Voronoi cell contains the point \p p.
   */
  Vertex_handle_rt find_seed_near_point(const double* p) {
    // In order to be compatible with the symbolic
    // perturbation, if the nearest neighbor is
    // non-unique, we need to return the one of
    // lowest index (because in case of several seeds
    // at equal distance, the one of lowest index
    // is guaranteed to have the facet in its Voronoi
    // cell from the point of view of symbolic
    // perturbation).
    Point pt(p[0], p[1], p[2]);
    return rt_->nearest_power_vertex(pt);
  }

  /********************************************************************/
 public:
  /**
   * \brief Gets the input mesh.
   */
  const GEO::Mesh* mesh() const { return mesh_; }

  /**
   * \brief Gets the input mesh.
   */
  GEO::Mesh* mesh() { return mesh_; }

  /**
   * \brief Gets the Delaunay triangulation.
   */
  RegularTriangulationNN* delaunay() { return rt_; }

  /**
   * \brief Sets the Delaunay triangulation.
   */
  void set_delaunay(RegularTriangulationNN* rt) { rt_ = rt; }

  /**
   * \brief Sets the input mesh.
   */
  void set_mesh(GEO::Mesh* mesh) { mesh_ = mesh; }
  /**
   * \brief Sets the facets range.
   * \details Computations can be restricted to a contiguous facet range.
   * \param[in] facets_begin first facet in the range.
   * \param[in] facets_end one position past the last facet in the range.
   */
  void set_facets_range(index_t facets_begin, index_t facets_end) {
    geo_debug_assert(facets_end >= facets_begin);
    facets_begin_ = facets_begin;
    facets_end_ = facets_end;
  }

  /**
   * \brief Gets the current cell.
   * \details The current cell corresponds to the
   *  intersection between the current tetrahedron
   *  and the Voronoi cell of the current seed.
   *  Can be used in volumetric traversals (and not in
   *  surfacic traversals).
   */
  const Polyhedron& current_polyhedron() const { return *current_polyhedron_; }

  /**
   * \brief Sets symbolic mode.
   * \details If exact mode is active, symbolic mode is enforced.
   * \param[in] x if set, the symbolic representation of the intersections
   *  are computed.
   */
  void set_symbolic(bool x) {
    symbolic_ = x;
    // exact mode requires symbolic mode.
    if (exact_) {
      symbolic_ = true;
    }
  }
  /**
   * \brief Tests whether symbolic mode is active.
   */
  bool symbolic() const { return symbolic_; }
  /**
   * \brief Specifies whether exact predicates should be used.
   * \details If exact predicates are used, symbolic mode is ensured.
   * \param[in] x if set, exact predicates are used.
   */
  void set_exact_predicates(bool x) {
    exact_ = x;
    // exact mode requires symbolic mode.
    if (exact_) {
      symbolic_ = true;
    }
  }
  /**
   * \brief Tests whether exact predicates are used.
   */
  bool exact_predicates() const { return exact_; }

  /**
   * \brief Specifies whether radius of security should be enforced.
   */
  void set_check_SR(bool x) { check_SR_ = x; }

  /**
   * \brief Tests whether radius of security is enforced.
   * \retval true if radius of security test is used.
   * \retval false otherwise.
   */
  bool check_SR() const { return check_SR_; }

  /**
   * \brief Gets the PointAllocator.
   * \return a pointer to the PointAllocator, used
   *  to create the new vertices generated by
   *  intersections.
   */
  PointAllocator* point_allocator() { return &intersections_; }

 protected:
  /**
   * \name Adapter classes for surfacic computation
   * @{
   */

  /**
   * \brief Adapter class used internally to implement for_each_polygon()
   * \details Overrides constness checks, to allow using temporaries as
   *   argument of for_each_xxx().
   * \tparam ACTION the user action class.
   */
  template <class ACTION>
  class PolygonAction {
   public:
    /**
     * \brief Creates a new PolygonAction around a user ACTION instance.
     * \param[in] do_it the user ACTION instance
     */
    PolygonAction(const ACTION& do_it) : do_it_(do_it) {}

    /**
     * \brief Callback called for each polygon.
     * \details Routes the callback to the wrapped user action class.
     * \param[in] v index of current Delaunay seed
     * \param[in] f index of current mesh facet
     * \param[in] P intersection between current mesh facet
     *  and the Voronoi cell of \p v
     */
    void operator()(const index_t v, const index_t f, const Polygon& P) const {
      // std::cout << "In PolygonAction operator..." << std::endl;
      GEO::geo_argused(f);
      // std::cout << "start casting action" << std::endl;
      const_cast<ACTION&>(do_it_)(v, P);
      // std::cout << "finish casting action" << std::endl;
    }

   protected:
    const ACTION& do_it_;
  };

  /**
   * \brief Adapter class used internally to implement
   *  for_each_triangle().
   * \details Overrides constness checks, to allow using temporaries as
   * argument of for_each_xxx().
   * \tparam ACTION the user action class
   */
  template <class ACTION>
  class TriangleAction {
   public:
    /**
     * \brief Creates a new TriangleAction that wraps a
     *  user ACTION instance.
     * \param[in] do_it the user ACTION instance
     */
    TriangleAction(const ACTION& do_it) : do_it_(do_it) {}

    /**
     * \brief Callback called for each integration simplex.
     * \details Decomposes the polygon \p P into triangles and
     *  calls the callback of the wrapped user action class
     *  for each triangle.
     * \param[in] v index of current Delaunay seed
     * \param[in] f index of current mesh facet
     * \param[in] P intersection between current mesh facet and
     *  the Voronoi cell of \p v
     */
    void operator()(index_t v, index_t f, const Polygon& P) const {
      GEO::geo_argused(f);
      // if (v == 457) {
      //     for(index_t i = 0; i < P.nb_vertices(); i++) {
      //         const double* ve = P.vertex(i).point();
      //         std::cout << "v: " << v << " has polygon vertex (" << ve[0] <<
      //         "," << ve[1] << "," << ve[2] << ")" << std::endl;
      //     }
      // }

      for (index_t i = 1; i + 1 < P.nb_vertices(); i++) {
        // // handle degenerated cases
        // if (is_same(P.vertex(0), P.vertex(i)))
        //     continue;

        const_cast<ACTION&>(do_it_)(v, P.vertex(0), P.vertex(i),
                                    P.vertex(i + 1));
      }
    }

   protected:
    static bool is_same(const double* p0, const double* p1) {
      if (std::abs(p0[0] - p1[0]) <= SCALAR_ZERO &&
          std::abs(p0[1] - p1[1]) <= SCALAR_ZERO &&
          std::abs(p0[2] - p1[2]) <= SCALAR_ZERO) {
        return true;
      } else {
        return false;
      }
    }

   protected:
    const ACTION& do_it_;
  };

  /**
   * @}
   * \name Adapter classes for volumetric computation
   * @{
   */

  /**
   * \brief Adapter class used internally to implement
   *  for_each_polyhedron()
   * \details Overrides constness checks, to allow using temporaries as
   *  argument of for_each_xxx()
   * \tparam ACTION the user action class
   */
  template <class ACTION>
  class PolyhedronAction {
   public:
    /**
     * \brief Creates a new PolyhedronAction that wraps
     *  a user ACTION instance.
     * \param[in] do_it the user ACTION instance
     */
    PolyhedronAction(const ACTION& do_it) : do_it_(do_it) {}

    /**
     * \brief Callback called for each polyhedron
     * \details Routes the callback to the wrapped user action class.
     * \param[in] v index of current Delaunay seed
     * \param[in] t index of current mesh tetrahedron
     * \param[in] C intersection between current mesh tetrahedron
     *  and the Voronoi cell of \p v
     */
    void operator()(index_t v, index_t t, const Polyhedron& C) const {
      const_cast<ACTION&>(do_it_)(v, t, C);
    }

   protected:
    const ACTION& do_it_;
  };

  /**
   * \brief Adapter class used internally to implement
   *  for_each_volumetric_integration_simplex()
   * \details Overrides constness checks, to allow using temporaries as
   *  argument of for_each_xxx().
   * \tparam ACTION the user action class. It needs to implement:
   *  operator()(index_t v, signed_index_t v_adj,
   *    index_t t, signed_index_t t_adj,
   *    const Vertex& v1, const Vertex& v2, const Vertex& v3
   *  )
   *  where the parameters are as follows:
   *    - v is the index of the current Voronoi cell
   *    (or Delaunay vertex)
   *    - v_adj is the index of the Voronoi cell adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    adjacent to v or -1 if current face is a tetrahedron facet
   *    - t is the index of the current tetrahedron
   *    - t_adj is the index of the tetrahedron adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    - v1,v2 and v3 are the three vertices of the facet on the
   *    border of the restricted Voronoi cell.
   */
  template <class ACTION>
  class VolumetricIntegrationSimplexAction {
   public:
    /**
     * \brief Creates a new VolumetricIntegrationSimplexAction
     *  that wraps a user ACTION instance.
     * \param[in] do_it the user ACTION instance
     * \param[in] visit_inner_tets if set, all the tetrahedron-cell
     *  intersections are visited, else only tetrahedra on the border
     *  of the restricted Voronoi cell are visited. Since all the
     *  visited triangles are connected to the current Voronoi seed
     *  by a tetrahedron, the computed volume is the same
     *  in both cases.
     * \param[in] coherent_triangles if set, this ensures that
     *  the polygonal facets of the cells are always triangulated
     *  in a coherent manner when seen from two different cells.
     *  For instance, it is required if a tetrahedral mesh is
     *  reconstructed.
     */
    VolumetricIntegrationSimplexAction(const ACTION& do_it,
                                       bool visit_inner_tets = false,
                                       bool coherent_triangles = false)
        : do_it_(do_it),
          visit_inner_tets_(visit_inner_tets),
          coherent_triangles_(coherent_triangles) {}

    /**
     * \brief Callback called for each polyhedron
     * \details Routes the callback to the wrapped user action class.
     * \param[in] v index of current Delaunay seed
     * \param[in] t index of current mesh tetrahedron
     * \param[in] C intersection between current mesh tetrahedron
     *  and the Voronoi cell of \p v
     */
    void operator()(index_t v, index_t t, const Polyhedron& C) const {
      for (index_t cv = 0; cv < C.max_v(); ++cv) {
        signed_index_t ct = C.vertex_triangle(cv);
        if (ct == -1) {
          continue;
        }
        geo_debug_assert(C.triangle_is_used(index_t(ct)));

        signed_index_t adjacent = C.vertex_id(cv);
        signed_index_t v_adj = -1;
        signed_index_t t_adj = -1;

        if (adjacent < 0) {
          // Negative adjacent indices correspond to
          // tet-tet links (ignored when we want to triangulate
          // the border of the restricted Voronoi cell while
          // ignoring internal structures).
          if (!visit_inner_tets_) {
            continue;
          }
          t_adj = -adjacent - 1;
        } else if (adjacent > 0) {
          // Positive adjacent indices correspond to
          // Voronoi seed - Voronoi seed link
          v_adj = adjacent - 1;
        }
        // and adjacent indicex equal to zero corresponds
        // to tet on border.

        Polyhedron::Corner c1(index_t(ct),
                              index_t(C.find_triangle_vertex(index_t(ct), cv)));

        // If required, ensure that two polygonal facets
        // seen from two different volumetric cells will
        // be triangulated coherently.
        if (coherent_triangles_) {
          move_to_first_corner_of_facet(C, c1, v);
        }

        const Vertex& v1 = C.triangle_dual(c1.t);

        Polyhedron::Corner c2 = c1;
        C.move_to_next_around_vertex(c2);
        geo_debug_assert(c2 != c1);

        Polyhedron::Corner c3 = c2;
        C.move_to_next_around_vertex(c3);
        geo_debug_assert(c3 != c1);
        do {
          const Vertex& v2 = C.triangle_dual(c2.t);
          const Vertex& v3 = C.triangle_dual(c3.t);
          const_cast<ACTION&>(do_it_)(v, v_adj, t, t_adj, v1, v2, v3);
          c2 = c3;
          C.move_to_next_around_vertex(c3);
        } while (c3 != c1);
      }
    }

    /**
     * \brief Finds the first corner of a facet in a Polyhedron.
     * \details This function is used to ensure that a facet is
     *  triangulated coherently when seen from two different
     *  volumetric cells, by generating a fan of triangles
     *  that radiates from the first corner. The global order
     *  used to find the first
     *  corner is defined by the function symbolic_compare().
     *
     * \param[in] C the Polyhedron
     * \param[in,out] c a corner of the facet, replaced by the
     *  first corner of the facet on exit.
     * \param[in] center_vertex_id index of the current Voronoi seed
     *  (needed to determine the full symbolic information in the
     *  vertices).
     */
    void move_to_first_corner_of_facet(const Polyhedron& C,
                                       Polyhedron::Corner& c,
                                       index_t center_vertex_id) const {
      Polyhedron::Corner first = c;
      Polyhedron::Corner cur = c;
      do {
        if (symbolic_compare(C.triangle_dual(cur.t), C.triangle_dual(c.t),
                             center_vertex_id)) {
          c = cur;
        }
        C.move_to_next_around_vertex(cur);
      } while (cur != first);
    }

    /**
     * \brief Compares the symbolic information of two vertices
     *  in such a way that a global order is defined.
     * \details This function is used to ensure that a facet is
     *  triangulated coherently when seen from two different
     *  volumetric cells (it uniquely determines the "first" vertex).
     * \param[in] p1 first vertex to compare
     * \param[in] p2 second vertex to compare
     * \param[in] center_vertex_id index of the current Voronoi seed
     *  (needed to determine the full symbolic information in
     *   \p p1 and \p p2).
     * \return true if p1 is before p2 in the global order,
     *  false otherwise.
     */
    static bool symbolic_compare(const Vertex& p1, const Vertex& p2,
                                 index_t center_vertex_id) {
      GEO::signed_quadindex K1(signed_index_t(center_vertex_id), p1.sym()[0],
                               p1.sym()[1], p1.sym()[2]);
      GEO::signed_quadindex K2(signed_index_t(center_vertex_id), p2.sym()[0],
                               p2.sym()[1], p2.sym()[2]);
      return K1 < K2;
    }

   protected:
    const ACTION& do_it_;
    bool visit_inner_tets_;
    bool coherent_triangles_;
  };

  /**
   * \brief Adapter class used internally to implement
   *  for_each_tetrahedron()
   * \details Overrides constness checks, to allow using temporaries as
   *  argument of for_each_xxx()
   * \tparam ACTION the user action class. It needs to implement:
   *  operator()(index_t v, signed_index_t v_adj,
   *    index_t t, index_t t_adj,
   *    const Vertex& v0, const Vertex& v1,
   *    const Vertex& v2, const Vertex& v3
   *  )
   *  where the parameters are as follows:
   *    - v is the index of the current Voronoi cell
   *    (or Delaunay vertex)
   *    - v_adj is the index of the Voronoi cell adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    adjacent to v or -1 if current face is a tetrahedron facet
   *    - t is the index of the current tetrahedron
   *    - t_adj is the index of the tetrahedron adjacent to t accros
   *    facet (\p v1, \p v2, \p v3) or -1 if it does not exists
   *    - v0,v1,v2 and v3 are the four vertices of tetrahedron.
   */
  template <class ACTION>
  class TetrahedronAction {
   public:
    /**
     * \brief Creates a new TetrahedronAction that wraps
     *  a user ACTION instance.
     * \param[in] do_it the user ACTION instance
     */
    TetrahedronAction(const ACTION& do_it) : do_it_(do_it) {}

    /**
     * \brief Callback called for each polyhedron
     * \details Routes the callback to the wrapped user action class.
     * \param[in] v index of current Delaunay seed
     * \param[in] t index of current mesh tetrahedron
     * \param[in] C intersection between current mesh tetrahedron
     *  and the Voronoi cell of \p v
     */
    void operator()(index_t v, index_t t, const Polyhedron& C) const {
      // Find a vertex of the current cell,
      // that will be used as the 'origin'
      // vertex
      const Vertex* v0 = nullptr;
      index_t t0;
      for (t0 = 0; t0 < C.max_t(); ++t0) {
        if (C.triangle_is_used(t0)) {
          v0 = &C.triangle_dual(t0);
          break;
        }
      }

      // If current cell is empty, return
      if (v0 == nullptr) {
        return;
      }

      for (index_t cv = 0; cv < C.max_v(); ++cv) {
        signed_index_t ct = C.vertex_triangle(cv);
        if (ct == -1) {
          continue;
        }
        geo_debug_assert(C.triangle_is_used(index_t(ct)));

        signed_index_t adjacent = C.vertex_id(cv);
        signed_index_t v_adj = -1;
        signed_index_t t_adj = -1;

        if (adjacent < 0) {
          // Negative adjacent indices correspond to
          // tet-tet links
          t_adj = -adjacent - 1;
        } else if (adjacent > 0) {
          // Positive adjacent indices correspond to
          // Voronoi seed - Voroni seed link
          v_adj = adjacent - 1;
        }
        // and adjacent indicex equal to zero corresponds
        // to tet on border.

        Polyhedron::Corner c1(index_t(ct),
                              C.find_triangle_vertex(index_t(ct), cv));

        // If the current facet is incident to
        // the origin vertex, then skip it (else
        // it would generate flat tetrahedra)
        if (facet_is_incident_to_vertex(C, c1, t0)) {
          continue;
        }

        const Vertex& v1 = C.triangle_dual(c1.t);

        Polyhedron::Corner c2 = c1;
        C.move_to_next_around_vertex(c2);
        geo_debug_assert(c2 != c1);

        Polyhedron::Corner c3 = c2;
        C.move_to_next_around_vertex(c3);
        geo_debug_assert(c3 != c1);
        do {
          const Vertex& v2 = C.triangle_dual(c2.t);
          const Vertex& v3 = C.triangle_dual(c3.t);
          const_cast<ACTION&>(do_it_)(v, v_adj, t, t_adj, *v0, v1, v2, v3);
          c2 = c3;
          C.move_to_next_around_vertex(c3);
        } while (c3 != c1);
      }
    }

   protected:
    /**
     * \brief Tests whether a Polyhedron facet is incident
     *  to a vertex.
     * \param[in] C the Polyhedron
     * \param[in] c a corner of the facet
     * \param[in] t the index of the vertex in dual form (in other
     *  words, a triangle index).
     * \return true if the facet incident to corner \p c
     *  is also incident to the vertex dual to \p t, false otherwise
     */
    bool facet_is_incident_to_vertex(const Polyhedron& C, Polyhedron::Corner& c,
                                     index_t t) const {
      Polyhedron::Corner first = c;
      Polyhedron::Corner cur = c;
      do {
        if (cur.t == t) {
          return true;
        }
        C.move_to_next_around_vertex(cur);
      } while (cur != first);
      return false;
    }

   protected:
    const ACTION& do_it_;
  };

 protected:
  RegularTriangulationNN* rt_;
  GEO::Mesh* mesh_;

  PointAllocator intersections_;
  Polygon* current_polygon_;
  Polygon P1, P2;
  std::vector<Vertex_handle_rt> neighbors_;
  index_t current_facet_;
  index_t current_seed_;
  Vertex_handle_rt current_seed_handle_;
  Polyhedron* current_polyhedron_;
  index_t current_tet_;

  // // For optimized get_neighbors().
  // GEO::signed_index_t cur_stamp_;
  // GEO::vector<signed_index_t> stamp_;

  bool symbolic_;
  bool check_SR_;
  bool exact_;

  // though we have weight, dimension should still be 3
  coord_index_t dimension_;

  static const index_t UNSPECIFIED_RANGE = index_t(-1);

  index_t facets_begin_;
  index_t facets_end_;

  index_t tets_begin_;
  index_t tets_end_;

  bool connected_components_priority_;

 private:
  /**
   * \brief Forbids construction from copy.
   */
  GenRestrictedPowerDiagram(const thisclass& rhs);

  /**
   * \brief Forbids assignment.
   */
  thisclass& operator=(const thisclass& rhs);
};

// namespace matfp {
//     /**
//      * \brief Symbolic representation of a RestrictedVoronoiDiagram vertex.
//      */
//     typedef matfpGen::SymbolicVertex SymbolicVertex;
// }

}  // namespace matfp
