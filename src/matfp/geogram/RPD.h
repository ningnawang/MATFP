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

#include <geogram/basic/attributes.h>
#include <geogram/basic/common.h>
#include <geogram/basic/counted.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/smart_pointer.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh.h>

#include "matfp/Types/CommonTypes.h"
#include "matfp/Types/RegularTriangulation.h"
#include "matfp/geogram/RPD_callback.h"
#include "matfp/geogram/RPD_mesh_builder.h"
#include "matfp/geogram/generic_RPD.h"

namespace matfp {
using namespace GEO;

using GEOGen::PointAllocator;

class GEOGRAM_API RestrictedPowerDiagram : public GEO::Counted {
 public:
  static RestrictedPowerDiagram* create(RegularTriangulationNN* rt,
                                        GEO::Mesh* mesh,
                                        const double* R3_embedding,
                                        index_t R3_embedding_stride);

  static RestrictedPowerDiagram* create(RegularTriangulationNN* rt,
                                        GEO::Mesh* mesh) {
    return create(
        rt, mesh,
        (mesh->vertices.nb() > 0) ? mesh->vertices.point_ptr(0) : nullptr,
        mesh->vertices.dimension());
  }

  /**
   * \brief Computes the restricted Power diagram and stores it
   *  in a mesh.
   * \param[out] M the computed restricted Voronoi diagram
   * \param[in] dim if different from 0, use only the
   *  first dim coordinates
   * \param[in] cell_borders_only in volumetric mode, computes only
   *  the surfacic borders of the volumetric cells (for visualization
   *  purpose)
   * \param[in] integration_simplices in volumetric mode, if set,
   *  the generated tetrahedra systematically have the Voronoi seed
   *  as the first vertex. As a consequence, the mesh is not necessarily
   *  geometrically correct (it may have inverted elements), but it is
   *  algebraically correct (the sum of signed volumes corresponds the
   *  the total volume of each cell).
   */
  virtual void compute_RPD(
      GEO::Mesh& M,
      std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_seed_adj,
      std::map<GEO::index_t, std::set<GEO::index_t>>* facet_adj_seeds,
      GEO::coord_index_t dim = 0, bool cell_borders_only = false,
      bool integration_simplices = false, bool is_parallel = true) = 0;

  /**
   * \brief Gets the dimension used by this RestrictedPowerDiagram.
   */
  coord_index_t dimension() const { return dimension_; }

  /**
   * \brief Gets the Delaunay triangulation.
   */
  RegularTriangulationNN* delaunay() { return rt_; }

  /**
   * \brief Sets the Regular Triangulation.
   */
  virtual void set_delaunay(RegularTriangulationNN* rt);

  /**
   * \brief Tests whether volumetric mode is used.
   */
  bool volumetric() const { return volumetric_; }

  /**
   * \brief Sets volumetric mode.
   * \param[in] x if true, volumetric mode is used, otherwise
   *  surfacic mode is used.
   */
  virtual void set_volumetric(bool x) = 0;

  /**
   * \brief Invokes a user callback for each intersection polygon
   *  of the restricted Voronoi diagram (surfacic mode only).
   * \details Each intersection polygon is defined as the intersection
   *  between a Voronoi cell and a triangle.
   * \param[in] callback the set of user callbacks, as an instance of a
   *  class derived from RPDPolygonCGALCallback.
   * \param[in] symbolic if true, generate symbolic information in the
   *  vertices
   * \param[in] connected_comp_priority if true, generate polyhedron
   *  intersections associated with the same Voronoi seed in order.
   * \param[in] parallel if true, tentatively parallelize computation.
   */
  virtual void for_each_polygon(matfp::RPDPolygonCGALCallback& callback,
                                bool symbolic = true,
                                bool connected_comp_priority = true,
                                bool parallel = false) = 0;

  /**
   * \brief Invokes a user callback for each intersection polyhedron
   *  of the restricted Voronoi diagram (volumetric mode only).
   * \details Each intersection polyhedron is defined as the intersection
   *  between a Voronoi cell and a tetrahedron.
   * \param[in] callback the set of user callbacks, as an instance of a
   *  class derived from RVDPolyhedronCallback.
   * \param[in] symbolic if true, generate symbolic information in the
   *  vertices
   * \param[in] connected_comp_priority if true, generate polyhedron
   *  intersections associated with the same Voronoi seed in order.
   * \param[in] parallel if true, tentatively parallelize computation.
   */
  virtual void for_each_polyhedron(matfp::RPDPolyhedronCGALCallback& callback,
                                   bool symbolic = true,
                                   bool connected_comp_priority = true,
                                   bool parallel = false) = 0;

  /**
   * \brief Specifies whether the "radius of security"
   *  criterion should be enforced.
   */
  virtual void set_check_SR(bool x) = 0;
  /**
   * \brief Specifies whether exact predicates should
   *  be used.
   */
  virtual void set_exact_predicates(bool x) = 0;
  /**
   * \brief Partitions the mesh and creates
   *  local storage for multithreaded implementation.
   */
  virtual void create_threads() = 0;

  /**
   * \brief Deletes all the local storage associated
   *  with the threads.
   */
  virtual void delete_threads() = 0;
  /**
   * \brief Restricts surfacic computations to a part of the input mesh.
   * \details The part of the input mesh should be specified as
   *    a contiguous range of facet indices.
   * \param[in] facets_begin first facet in the range
   * \param[in] facets_end one past last facet in the range
   */
  virtual void set_facets_range(index_t facets_begin, index_t facets_end) = 0;
  /**
   * \brief Gets the PointAllocator.
   * \return a pointer to the PointAllocator, used
   *  to create the new vertices generated by
   *  intersections.
   */
  virtual GEOGen::PointAllocator* point_allocator() = 0;

 protected:
  /**
   * \brief This constructor is never called directly.
   * \details Use one of the three versions of create() instead.
   */
  RestrictedPowerDiagram(RegularTriangulationNN* rt, Mesh* mesh,
                         const double* R3_embedding,
                         index_t R3_embedding_stride);

  /**
   * \brief RestrictedPowerDiagram destructor
   */
  virtual ~RestrictedPowerDiagram();

 protected:
  coord_index_t dimension_;
  RegularTriangulationNN* rt_;
  GEO::Mesh* mesh_;
  const double* R3_embedding_base_;
  index_t R3_embedding_stride_;
  signed_index_t facets_begin_;
  signed_index_t facets_end_;
  signed_index_t tets_begin_;
  signed_index_t tets_end_;
  bool volumetric_;
};

/** \brief Smart pointer to a RestrictedPowerDiagram object */
typedef GEO::SmartPointer<RestrictedPowerDiagram> RestrictedPowerDiagram_var;

}  // namespace matfp