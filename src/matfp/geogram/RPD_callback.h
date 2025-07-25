/*
 *  Copyright (c) 2010-2017, ALICE project, Inria
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

#include <geogram/basic/common.h>
#include <geogram/voronoi/RVD_callback.h>

#include "matfp/geogram/RPD_mesh_builder.h"
#include "matfp/geogram/generic_RPD_cell.h"
#include "matfp/geogram/generic_RPD_polygon.h"
#include "matfp/geogram/generic_RPD_vertex.h"  // for GEO::RVDCallback
// #include <geogram/voronoi/generic_RVD_vertex.h> // for GEO::RVDCallback
#include <geogram/basic/numeric.h>
#include <geogram/mesh/mesh.h>

// namespace GEO {
// namespace Process {
// class SpinLockArray;
// }
// }  // namespace GEO

/**
 * \file geogram/voronoi/RVD_callback.h
 * \brief Declaration of base types for implementing user-defined
 *  code that queries the cells of a restricted Voronoi diagram.
 */

namespace matfp {
using namespace GEO;

/**
 * \brief Baseclass for user functions called for each
 *  polygon of a surfacic restricted Voronoi diagram.
 * \details A surfacic restricted Voronoi diagram is the
 *  intersection between Voronoi cells and the triangles of
 *  a given triangulated mesh. The member functions of this
 *  class are called for each intersection between a Voronoi
 *  cell and a triangle.
 */

class GEOGRAM_API RPDPolygonCGALCallback : public GEO::RVDCallback {
 public:
  /**
   * \brief PolyhedronCallback constructor.
   */
  RPDPolygonCGALCallback();

  /**
   * \brief PolyhedronCallback destructor.
   */
  virtual ~RPDPolygonCGALCallback();

  /**
   * \copydoc RVDCallback::begin()
   */
  virtual void begin();

  /**
   * \copydoc RVDCallback::end()
   */
  virtual void end();

  /**
   * \brief The default callback called for each polygon
   * \param[in] v index of current Delaunay seed
   * \param[in] t index of current mesh triangle
   * \param[in] C intersection between current mesh triangle
   *  and the Voronoi cell of \p v
   */
  virtual void operator()(GEO::index_t v, GEO::index_t t,
                          const matfp::PolygonCGAL& C) const;
};

/***************************************************************/

/**
 * \brief Baseclass for user functions called for each
 *  polyhedron of a volumetric restricted Voronoi diagram.
 * \details A volumetric restricted Voronoi diagram is the
 *  intersection between Voronoi cells and the tetrahedra of
 *  a given tetrahedral mesh. The member functions of this
 *  class are called for each intersection between a Voronoi
 *  cell and a tetrahedron.
 */
class GEOGRAM_API RPDPolyhedronCGALCallback : public RVDCallback {
 public:
  /**
   * \brief PolyhedronCallback constructor.
   */
  RPDPolyhedronCGALCallback();

  /**
   * \brief PolyhedronCallback destructor.
   */
  virtual ~RPDPolyhedronCGALCallback();

  /**
   * \copydoc RVDCallback::begin()
   */
  virtual void begin();

  /**
   * \copydoc RVDCallback::end()
   */
  virtual void end();

  /**
   * \brief The default callback called for each polyhedron
   * \details This default implementation routes the callback to the
   *  begin_polyhedron_internal(), end_polyhedron_internal(),
   *  begin_facet_internal(), end_facet_internal() and vertex_internal()
   *  functions (that in turn route the callbacks to their without
   *  "_internal" counterparts).
   * \param[in] v index of current Delaunay seed
   * \param[in] t index of current mesh tetrahedron
   * \param[in] C intersection between current mesh tetrahedron
   *  and the Voronoi cell of \p v
   */
  virtual void operator()(GEO::index_t v, GEO::index_t t,
                          const matfp::ConvexCellCGAL& C) const;

  /**
   * \brief Called at the beginning of each intersection polyhedron.
   * \details Each intersection polyhedron is defined as the intersection
   *   between a Voronoi cell and a tetrahedron.
   * \param[in] seed index of the seed associated with the Voronoi cell
   * \param[in] tetrahedron index of the tetrahedron
   */
  virtual void begin_polyhedron(GEO::index_t seed, GEO::index_t tetrahedron);

  /**
   * \brief Called at the beginning of each facet of each intersection
   *   polyhedron.
   * \details A facet can be a subset of either a bisector (defined by
   *   two seeds), or it can be a subset of a facet of a tetrahedron.
   * \param[in] facet_seed if the facet corresponds to a bisector,
   *   the index of the seed that defines the bisector, or index_t(-1)
   *   otherwise
   * \param[in] facet_tet if the facet corresponds to a facet
   *   of the current tetrahedron, the index of the tetrahedron
   *   adjacent to that facet, or index_t(-1) otherwise
   */
  virtual void begin_facet(GEO::index_t facet_seed, GEO::index_t facet_tet);

  /**
   * \brief Called for each vertex of the current facet.
   * \param[in] geometry a pointer to the coordinates of the vertex
   * \param[in] symb the symbolic representation of the vertex
   */
  virtual void vertex(const double* geometry,
                      const matfp::SymbolicVertex& symb);

  /**
   * \brief Called at the end of each polyhedron facet.
   */
  virtual void end_facet();

  /**
   * \brief Called at the end of each polyhedron.
   */
  virtual void end_polyhedron();

  /**
   * \brief Gets the index of the tetrahedron that corresponds to the
   *  current polyhedron.
   * \return The index of the tetrahedron that corresponds to the
   *  current polyhedron.
   * \details The current polyhedron is the intersection between a Voronoi
   *  cell and a tetrahedron.
   */
  GEO::index_t tet() const { return simplex(); }

  /**
   * \brief Gets the index of the seed that defines the bisector on which
   *  the current facet lies, or index_t(-1).
   * \return The index of the seed that defines the bisector on which
   *  the current facet lies, or index_t(-1).
   * \details Each facet is either on a bisector or on a tetrahedron
   *  facet. If the current facet is on a bisector, it is defined by
   *  seed() and facet_seed(), otherwise facet_seed() returns index_t(-1).
   */
  GEO::index_t facet_seed() const { return facet_seed_; }

  /**
   * \brief Gets the index of the tetrahedron adjacent to the current
   *  facet or index_t(-1) if there is no such facet.
   * \return the index of the tetrahedron adjacent to the current
   *  facet or index_t(-1).
   * \details Each facet is either on a bisector or on a tetrahedron
   *  facet. If the current facet is on a tetrahedron facet, then it
   *  is defined by tet() and facet_tet(), otherwise
   *  facet_tet() returns index_t(-1).
   */
  GEO::index_t facet_tet() const { return facet_tet_; }

  /**
   * \brief Specifies whether internal tetrahedron facets should be
   *  removed.
   * \details If set, a single polyhedron is generated for each
   *  (connected component) of the restricted Voronoi cells. If not
   *  set (default), each tetrahedron-Voronoi cell intersection
   *  generates a new polyhedron.
   * \param[in] x true if internal facets should be removed, false
   *  otherwise
   */
  void set_simplify_internal_tet_facets(bool x) {
    simplify_internal_tet_facets_ = x;
  }

  /**
   * \brief Specifies whether Voronoi facets should be simplified.
   * \details By default, the computed Voronoi facets are composed
   *  of multiple polyhedra that correspond to the intersection with
   *  the tetrahedra of the input volume mesh. They can be simplified
   *  and replaced by a single polygon. This implies simplifying the
   *  internal tetrahedron facets and using a mesh.
   * \param[in] x true if Voronoi facets should be simplified,
   *  false otherwise.
   */
  void set_simplify_voronoi_facets(bool x) {
    simplify_voronoi_facets_ = x;
    if (x) {
      set_simplify_internal_tet_facets(true);
      set_use_mesh(true);
    }
  }

  /**
   * \brief Specifies whether boundary facets should be simplified.
   * \details By default, the intersection between a Voronoi cell and
   *  the boundary is possibly composed of multiple polygons, that
   *  correspond to the initial polygons of the boundary. They can be
   *  simplified as a single polygon per Voronoi cell. This implies
   *  simplifying the internal tetrahedron facets, simplifying the
   *  Voronoi facets and using a mesh.
   * \param[in] x true if boundary facets should be simplified,
   *  false otherwise.
   * \param[in] angle_threshold an edge shared by two adjacent facets
   *  is suppressed if the angle between the facet normals is smaller
   *  than \p angle_threshold
   */
  void set_simplify_boundary_facets(bool x, double angle_threshold = 45.0) {
    simplify_boundary_facets_ = x;
    if (x) {
      set_simplify_voronoi_facets(true);
      simplify_boundary_facets_angle_threshold_ = angle_threshold;
    } else {
      simplify_boundary_facets_angle_threshold_ = 0.0;
    }
  }

  /**
   * \brief Specifies whether non-convex facets should be tessellated.
   * \param[in] x true if non-convex facets should be tessellated,
   *  false otherwise.
   * \details Only taken into account if set_use_mesh(true) was called.
   */
  void set_tessellate_non_convex_facets(bool x) {
    tessellate_non_convex_facets_ = x;
  }

  /**
   * \brief Specifies whether a mesh should be built for each
   *  traversed polyhedron.
   * \details The build mesh can then be modified (e.g., simplified)
   *  by overloading process_mesh().
   * \param[in] x true if a mesh should be build, false otherwise
   */
  void set_use_mesh(bool x);

  /**
   * \brief Sets the dimension of the internal mesh if need be.
   * \details This function is called automatically by
   *  RestrictedVoronoiDiagram::for_each_polyhedron().
   * \param[in] dim the dimension of the mesh (3 for 3d).
   */
  void set_dimension(GEO::index_t dim) {
    if (use_mesh_) {
      mesh_.vertices.set_dimension(dim);
    }
  }

 protected:
  /**
   * \brief Filters callbacks between operator() and client callbacks.
   * \details This is used to implement cells simplifications (remove
   *  internal boundaries and intersections with tetrahedra).
   * \see begin_polyhedron()
   */
  virtual void begin_polyhedron_internal(GEO::index_t seed,
                                         GEO::index_t tetrahedron);

  /**
   * \brief Filters callbacks between operator() and client callbacks.
   * \details This is used to implement cells simplifications (remove
   *  internal boundaries and intersections with tetrahedra).
   * \see begin_facet()
   */
  virtual void begin_facet_internal(GEO::index_t facet_seed,
                                    GEO::index_t facet_tet);

  /**
   * \brief Filters callbacks between operator() and client callbacks.
   * \details This is used to implement cells simplifications (remove
   *  internal boundaries and intersections with tetrahedra).
   * \see vertex()
   */
  virtual void vertex_internal(const double* geometry,
                               const matfp::SymbolicVertex& symb);

  /**
   * \brief Filters callbacks between operator() and client callbacks.
   * \details This is used to implement cells simplifications (remove
   *  internal boundaries and intersections with tetrahedra).
   * \see end_facet()
   */
  virtual void end_facet_internal();

  /**
   * \brief Filters callbacks between operator() and client callbacks.
   * \details This is used to implement cells simplifications (remove
   *  internal boundaries and intersections with tetrahedra).
   * \see end_polyhedron()
   */
  virtual void end_polyhedron_internal();

  /**
   * \brief If use_mesh is set, then this function is called for
   *  each generated mesh.
   * \details Default implementation simplifies the mesh based on
   *  sipmlify_xxx flags, then it calls user callbacks
   *  begin_facet(), end_facet(), vertex() for each facet of the mesh,
   *  as well as begin_polyhedron() and end_polyhedron() once per mesh.
   *  Derived classes may modify (e.g., simplify) the mesh before calling
   *  user callbacks.
   */
  virtual void process_polyhedron_mesh();

 protected:
  GEO::index_t facet_seed_;
  GEO::index_t facet_tet_;
  GEO::index_t last_seed_;

  bool simplify_internal_tet_facets_;
  bool simplify_voronoi_facets_;
  bool simplify_boundary_facets_;
  double simplify_boundary_facets_angle_threshold_;
  bool tessellate_non_convex_facets_;

  bool use_mesh_;
  bool facet_is_skipped_;

  GEO::Mesh mesh_;
  GEO::Attribute<matfp::SymbolicVertex> mesh_vertex_sym_;
  GEO::Attribute<GEO::index_t> mesh_facet_seed_;
  GEO::Attribute<GEO::index_t> mesh_facet_tet_;
  matfp::RPDVertexMap* vertex_map_;
  std::vector<GEO::index_t> base_current_facet_;
};

}  // namespace matfp