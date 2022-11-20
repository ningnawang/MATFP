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

#include <geogram/basic/assert.h>
#include <geogram/basic/attributes.h>
#include <geogram/basic/common.h>
#include <geogram/basic/process.h>
#include <geogram/delaunay/delaunay_nn.h>
#include <geogram/mesh/mesh.h>
#include <matfp/Logger.h>

// definition of GEOGen::small_set, GEOGen::ORIGINAL, GEOGen::PointAllocator
#include <geogram/voronoi/generic_RVD_vertex.h>

/**
 * \file geogram/voronoi/generic_RPD_vertex.h
 * \brief Types and utilities for manipulating vertices in geometric
 *  and symbolic forms in restricted Voronoi diagrams.
 * \note This file contains functions and classes used by the
 *  internal implementation of GEO::GenericVoronoiDiagram.
 *  Except some special uses, e.g. subclassing GEO::IntegrationSimplex,
 *  they are not meant to be used directly by client code.
 */

namespace matfp {
using namespace GEO;

using GEO::coord_index_t; /**< \brief type for coordinate indices      */
using GEO::Delaunay;      /**< \brief type for nD Delaunay triangulation */
using GEO::index_t;       /**< \brief type for indices (vertex and facet id) */
using GEO::Sign; /**< \brief type for signs (POSITIVE,ZERO or NEGATIVE) */
using GEO::signed_index_t; /**< \brief type for indices (can be <0)     */

using GEO::Mesh;

using GEOGen::PointAllocator;
using GEOGen::small_set;

/**
 * \brief A set of three integers that encodes the
 *  equation of a vertex in GenericVoronoiDiagram.
 *
 * \details
 * - Each positive entry i denotes the bisector of the segment that connects
 *   the center vertex to the i-th vertex (note that the center vertex
 *   needs to be stored elsewhere, but is known when a RPD is used,
 *   since we know which dual cell we are processing).
 *
 * - Each negative entry i denotes the i-th face in the boundary TriMesh.
 *   Note: indexing starts with 1 (resp. -1), 0 is kept for error codes.
 *
 * - There is some additional information for the following
 *   two configurations:
 *   - boundary vertex: (nb_boundary_facets = 3)
 *     the index of the boundary vertex is returned
 *     by get_boundary_vertex()
 *   - intersection between boundary edge and bisector:
 *     (nb_boundary_facets = 2)
 *     the indices v1,v2 of the extremities of the boundary edges
 *     are obtained by get_boundary_edge(v1,v2)
 *
 *  Doing so avoids recomputing vertices that we already know
 *  (and avoids numerical problems when the boundary surface has
 *  coplanar (or nearly coplanar) facets).
 *  It also allows using exact predicates (not implemented yet).
 *
 * \note This is an internal implementation class, not meant to be
 *  used by client code.
 */
class SymbolicVertex : public GEOGen::small_set<GEO::signed_index_t, 3> {
  /** \brief This class type */
  typedef SymbolicVertex thisclass;

  /** \brief The base class of this class */
  typedef GEOGen::small_set<GEO::signed_index_t, 3> baseclass;

 public:
  /**
   * \brief Creates an uninitialized SymbolicVertex.
   */
  SymbolicVertex() : v1_(0), v2_(0) {}

  /**
   * \brief Adds a bisector to the symbolic representation.
   */
  void add_bisector(index_t i) { baseclass::insert(signed_index_t(i) + 1); }

  /**
   * \brief Adds a boundary facet to the symbolic representation.
   */
  void add_boundary_facet(index_t i) {
    baseclass::insert(-signed_index_t(i) - 1);
  }

  /**
   * \brief Gets the number of boundary facets in the
   *  symbolic representation.
   */
  index_t nb_boundary_facets() const {
    index_t result = 0;
    for (auto it = baseclass::begin(); it != baseclass::end() && *it < 0;
         ++it) {
      result++;
    }
    return result;
  }

  /**
   * \brief Gets the number of bisectors in the symbolic representation.
   */
  index_t nb_bisectors() const {
    index_t result = 0;
    for (auto it = baseclass::end() - 1;
         it != baseclass::begin() - 1 && *it > 0; --it) {
      result++;
    }
    return result;
  }

  /**
   * \brief Casts a signed_index_t into an (unsigned) index_t.
   * \details In debug mode, throws an assertion failure
   *  exception whenever \p x is negative.
   */
  static index_t to_unsigned_int(signed_index_t x) {
    geo_debug_assert(x >= 0);
    return (index_t)(x);
  }

  /**
   * \brief Gets a bisector
   * \param[in] i local index of the bisector
   * \return the index of the Delaunay vertex that corresponds to
   *    the second extremity of the bisector
   * \pre i < nb_bisectors()
   */
  index_t bisector(signed_index_t i) const {
    geo_debug_assert(i < (signed_index_t)nb_bisectors());
    return to_unsigned_int((baseclass::end()[-1 - i]) - 1);
  }

  /**
   * \brief Gets a boundary facet
   * \param[in] i local index of the boundary facet
   * \return the index of the mesh facet
   * \pre i < nb_boundary_facets()
   */
  index_t boundary_facet(signed_index_t i) const {
    geo_debug_assert(i < (signed_index_t)nb_boundary_facets());
    return to_unsigned_int(-(baseclass::begin()[i]) - 1);
  }

  /**
   * \brief Tests whether a bisector is present in the
   *  symbolic representation this vertex.
   * \param[in] i global index of the bisector
   */
  bool has_bisector(index_t i) const {
    return baseclass::find(signed_index_t(i) + 1) != baseclass::end();
  }

  /**
   * \brief Tests whether a boundary facet is present in the
   *  symbolic representation of this vertex.
   * \param[in] i global index of the boundary facet
   */
  bool has_boundary_facet(index_t i) const {
    return baseclass::find(-signed_index_t(i) - 1) != baseclass::end();
  }

  /**
   * \brief Gets the global index of the boundary vertex that corresponds
   * to this vertex.
   * \pre nb_boundary_facets() == 3
   */
  index_t get_boundary_vertex() const {
    geo_debug_assert(nb_boundary_facets() == 3);
    geo_debug_assert(v1_ != 0);
    return v1_ - 1;
  }

  /**
   * \brief Gets the global indices of the boundary vertices that
   *  define the boundary edge on which this vertex is located.
   * \param[out] v1 index of the first extremity of the boundary edge
   * \param[out] v2 index of the second extremity of the boundary edge
   * \pre nb_boundary_facets() == 2
   */
  void get_boundary_edge(index_t& v1, index_t& v2) const {
    geo_debug_assert(nb_boundary_facets() == 2);
    geo_debug_assert(v1_ != 0);
    geo_debug_assert(v2_ != 0);
    v1 = v1_ - 1;
    v2 = v2_ - 1;
  }

  /**
   * \brief Sets the boundary vertex on which this vertex is located.
   * \param[in] v global index of the boundary vertex
   */
  void set_boundary_vertex(index_t v) {
    v1_ = v + 1;
    v2_ = 0;
  }

  /**
   * \brief Sets the boundary edge on which this vertex is located.
   * \param[in] v1 global index of the first boundary vertex
   * \param[in] v2 global index of the second boundary vertex
   */
  void set_boundary_edge(index_t v1, index_t v2) {
    v1_ = v1 + 1;
    v2_ = v2 + 1;
  }

  /**
   * \brief Copies a boundary edge from the symbolic representation
   *  of another vertex.
   */
  void copy_boundary_edge_from(const thisclass& rhs) {
    geo_debug_assert(rhs.nb_boundary_facets() == 2);
    geo_debug_assert(rhs.nb_bisectors() == 1);
    geo_debug_assert(rhs.v1_ > 0);
    geo_debug_assert(rhs.v2_ > 0);
    v1_ = rhs.v1_;
    v2_ = rhs.v2_;
  }

  /**
   * \brief Computes the symbolic representation of the intersection
   *  between a segment and a bisector.
   * \details Computes the intersection between
   *  the segment [\p v1, \p v2] and the bisector \p E
   *
   * \return false if there was a problem
   *  (happens sometimes in finite precision mode)
   */
  bool intersect_symbolic(const thisclass& v1, const thisclass& v2, index_t E) {
    // Compute the symbolic representation as the intersection
    // of three planes.
    sets_intersect(v1, v2, *this);
    // this computes the set of planes that contain
    // the edge [v1,v2]

    add_bisector(E);  // the intersection is on E.

    // Compute the symbolic representation as intersection between
    //    bisector and boundary edge
    // (it's redundant and less elegant than the representation
    //  as planes interactions,
    //  but we need this to handle degenerate configurations properly,
    //  and to use exact predicates with original boundary vertices
    //  coordinates).

    if (nb_boundary_facets() == 2) {
      // If *this is on the intersection of two boundary facets,
      // then *this is on
      // a boundary edge, and we need to retrieve the indices of the
      // two extremities of this boundary edge.

      index_t nb1 = v1.nb_boundary_facets();
      index_t nb2 = v2.nb_boundary_facets();
      if (nb1 == 3 && nb2 == 3) {
        // If v1 and v2 are boundary vertices,
        // then I is on the boundary
        // edge that connects v1 and v2
        set_boundary_edge(v1.get_boundary_vertex(), v2.get_boundary_vertex());
      } else if (nb1 == 2) {
        geo_debug_assert(nb_boundary_facets() == 2);
        // If v1 is on a boundary edge,
        // then I is on the same boundary edge as v1
        copy_boundary_edge_from(v1);
      } else if (nb2 == 2) {
        geo_debug_assert(nb_boundary_facets() == 2);
        // If v2 is on a boundary edge,
        // then I is on the same boundary edge as v2
        copy_boundary_edge_from(v2);
      }
    }

    // Sanity check: problem detected here, we
    // notify the caller that will use a workaround
    // (see clip_by_plane())
    if (baseclass::size() != 3) {
      return false;
    }
    return true;
  }

 private:
  index_t v1_;
  index_t v2_;
};

/**
 * \brief Internal representation of vertices
 *  in GenericVoronoiDiagram.
 * \details Vertex has both
 *  geometrical and symbolic representations.
 * \note This is an internal implementation class, not meant to be
 *  used by client code (except in some particular case, such as
 *  subclassing GEO::IntegrationSimplex).
 */
class Vertex {
  /** \brief This class type */
  typedef Vertex thisclass;

 public:
  /**
   * \brief Creates a new Vertex
   * \param[in] p geometric location at the vertex, shared with caller
   * \param[in] w weight
   * \param[in] f facet of the input mesh this Vertex comes from
   * \param[in] sym symbolic representation
   */
  Vertex(const double* p, double w, signed_index_t f, const SymbolicVertex& sym)
      : point_(p),
        weight_(w),
        f_(f),
        seed_(-1),
        sym_(sym),
        flags_(GEOGen::ORIGINAL) {}

  /**
   * \brief Creates a new Vertex
   * \param[in] p geometric location at the vertex, shared with caller
   * \param[in] w weight
   * \param[in] f facet of the input mesh this Vertex comes from
   */
  Vertex(const double* p, double w, signed_index_t f)
      : point_(p), weight_(w), f_(f), seed_(-1), flags_(GEOGen::ORIGINAL) {}

  /**
   * \brief Creates an uninitialized Vertex.
   */
  Vertex() : point_(nullptr), weight_(1.0), f_(-1), seed_(-1), flags_(0) {}

  /**
   * \brief Gets the geometric location at this Vertex.
   * \return a const pointer to the coordinates
   */
  const double* point() const { return point_; }

  /**
   * \brief Sets the geometric location at this vertex.
   * \param[in] p the geometric location, shared with caller
   */
  void set_point(const double* p) { point_ = p; }

  /**
   * \brief Gets Vertex weight.
   * \details Used by non-uniform centroidal
   *  Voronoi tesselation.
   */
  double weight() const { return weight_; }

  /**
   * \brief Sets the vertex weight.
   * \details Used by non-uniform centroidal
   * Voronoi tesselation..
   */
  void set_weight(double w) { weight_ = w; }

  /**
   * \brief Gets the adjacent seed.
   * \return the global index of the adjacent seed
   */
  signed_index_t adjacent_seed() const { return seed_; }

  /**
   * \brief Sets the adjacent seed.
   * \param[in] s the global index of the adjacent seed
   */
  void set_adjacent_seed(signed_index_t s) { seed_ = s; }

  /** Symbolic representation */

  /**
   * \brief Gets the symbolic representation.
   */
  const SymbolicVertex& sym() const { return sym_; }

  /**
   * \brief Gets the symbolic representation.
   */
  SymbolicVertex& sym() { return sym_; }

  /**
   * \brief Gets the adjacent facet.
   * \return the global index of the adjacent facet
   */
  signed_index_t adjacent_facet() const { return f_; }

  /**
   * \brief Sets the adjacent facet.
   * \param[in] f the global index of the adjacent facet
   */
  void set_adjacent_facet(signed_index_t f) { f_ = f; }

  /**
   * \brief Implicit conversion that accesses the geometric location.
   * \details With this implicit conversions, we can have template
   * arguments for RestrictedVoronoiDiagram that take
   * const double* as arguments instead of Vertices.
   * \return a const pointer to the coordinates
   */
  operator const double*() const { return point_; }

  /**
   * \brief Clears this Vertex.
   */
  void clear() {
    flags_ = 0;
    f_ = -1;
  }

  /**
   * \brief Sets an GEOGen::EdgeFlag in this Vertex.
   */
  void set_flag(GEOGen::EdgeFlag f) { flags_ |= f; }

  /**
   * \brief Resets an GEOGen::EdgeFlag in this Vertex.
   */
  void unset_flag(GEOGen::EdgeFlag f) { flags_ &= ~f; }

  /**
   * \brief Tests an GEOGen::EdgeFlag in this Vertex.
   */
  bool check_flag(GEOGen::EdgeFlag f) const { return (flags_ & f) != 0; }

  /**
   * \brief Copies adjacent facet and edge flags from another Vertex.
   */
  void copy_edge_from(const Vertex& rhs) {
    set_adjacent_facet(rhs.adjacent_facet());
    flags_ = rhs.flags_;
  }

  double dot_at(const double* p1, const double* p2, const double* ref,
                int DIM) {
    double r;
    for (coord_index_t c = 0; c < DIM; ++c) {
      r += (p1[c] - ref[c]) * (p2[c] - ref[c]);
    }
    return r;
  }

  double sq_dist(const double* p1, const double* p2, int DIM) {
    return dot_at(p1, p1, p2, DIM);
  }

  /**
   * \brief Computes the intersection between
   *    a segment and a bisector.
   * \details Computes the intersection between
   *  the segment [vq1, vq2] and the bisector
   *  of [p1,p2]..
   * \tparam DIM dimension, specified as a template
   *  argument for efficiency considerations
   */
  template <index_t DIM>
  void intersect_geom_new(GEOGen::PointAllocator& target_intersections,
                          const Vertex& vq1, const Vertex& vq2,
                          const double* p1, const double* p2) {
    const double l1 = sq_dist(p2, p1, DIM);
    const double R1 = p1[3] - p2[3];
    const double L1 = l1 + R1;

    // std::cout << "R1: " << R1 << std::endl;

    const double* q1 = vq1.point();
    const double* q2 = vq2.point();
    double* Ipoint = target_intersections.new_item();
    set_point(Ipoint);

    const double a10 = 2 * dot_at(p2, q1, p1, DIM);
    const double a11 = 2 * dot_at(p2, q2, p1, DIM);
    // const double Delta = ::fabs(a11 - a10);

    /*
     *       [ Lambda0 ]   [ -1 ]        [  a11 ]
     * Delta [         ] = [    ] * L1 + [      ]
     *       [ Lambda1 ]   [  1 ]        [ -a10 ]
     */
    double DeltaLambda1 = ::fabs(a11 - L1) * 0.5;
    double DeltaLambda2 = ::fabs(a10 - L1) * 0.5;
    const double Delta = DeltaLambda1 + DeltaLambda2;

    // logger().debug("DeltaLambda1: {}, DeltaLambda2: {}", DeltaLambda1,
    // DeltaLambda2);

    if (Delta > 1e-30) {
      DeltaLambda1 /= Delta;
      DeltaLambda2 /= Delta;
    } else {
      DeltaLambda1 = 0.5;
      DeltaLambda2 = 0.5;
    }
    for (coord_index_t c = 0; c < DIM; ++c) {
      Ipoint[c] = DeltaLambda1 * q1[c] + DeltaLambda2 * q2[c];
    }
    set_weight(DeltaLambda1 * vq1.weight() + DeltaLambda2 * vq2.weight());
  }

  /**
   * \brief Computes the intersection between
   *    a segment and a bisector.
   * \details Computes the intersection between
   *  the segment [vq1, vq2] and the bisector
   *  of [p1,p2]..
   * \tparam DIM dimension, specified as a template
   *  argument for efficiency considerations
   */
  template <index_t DIM>
  void intersect_geom(GEOGen::PointAllocator& target_intersections,
                      const Vertex& vq1, const Vertex& vq2, const double* p1,
                      const double* p2) {
    const double R1 = 0.5 * (p1[3] - p2[3]);
    const double* q1 = vq1.point();
    const double* q2 = vq2.point();
    double* Ipoint = target_intersections.new_item();
    set_point(Ipoint);
    double d = 0.0, l1 = 0.0, l2 = 0.0;
    for (coord_index_t c = 0; c < DIM; ++c) {
      double n = p1[c] - p2[c];
      d -= n * (p2[c] + p1[c]);
      l1 += q2[c] * n;
      l2 += q1[c] * n;
    }
    d = 0.5 * d;
    l1 = l1 + d + R1;
    l2 = l2 + d + R1;

    // l1 and l2 must be opposite sign
    // otherwise intersection is out of segment [q1,q2]
    if (l1 * l2 > 1e-30) {
      // logger().debug("!!!!!!!!!!!!!!! WRONG: l1: {}, l2: {}",  l1, l2);
      // log_and_throw("ERROR FOUND");
    }
    l1 = ::fabs(l1);
    l2 = ::fabs(l2);
    double l12 = l1 + l2;

    // l1 = ::fabs(l1 + d + R1);
    // l2 = ::fabs(l2 + d + R1);
    // double l12 = l1 + l2;

    // logger().debug("l1: {}, l2: {}", l1, l2);

    if (l12 > 1e-30) {
      l1 /= l12;
      l2 /= l12;
    } else {
      l1 = 0.5;
      l2 = 0.5;
    }
    for (coord_index_t c = 0; c < DIM; ++c) {
      Ipoint[c] = l1 * q1[c] + l2 * q2[c];
    }
    set_weight(l1 * vq1.weight() + l2 * vq2.weight());
  }

  /**
   * NOTE: only for RVD
   * \brief Computes the side of this vertex relative
   *  to a bisector.
   * \details This version is not exact.
   * \param[in] p1 first extremity of the bisector
   * \param[in] p2 second extremity of the bisector
   * \return POSITIVE if this vertex is on p1's side,
   *  NEGATIVE if this vertex is on p2's side, and ZERO
   *  if this vertex is on the bisector of [p1,p2].
   */
  template <index_t DIM>
  Sign side_fast(const double* p1, const double* p2) const {
    double r = 0.0;
    for (index_t c = 0; c < DIM; ++c) {
      r += GEO::geo_sqr(p2[c] - point()[c]);
      r -= GEO::geo_sqr(p1[c] - point()[c]);
    }
    return GEO::geo_sgn(r);
  }

 private:
  const double* point_;
  double weight_;

  /**
   * The facet adjacent to the edge
   * incident to this vertex.
   */
  signed_index_t f_;

  /**
   * indicates the seed of the bisector that generated the
   * edge that has this vertex and the previous one as
   * extremities (or -1 if border).
   */
  signed_index_t seed_;

  /** The symbolic representation of this vertex. */
  SymbolicVertex sym_;

  /**
   * Indicates the type of edge
   * (virtual, original or intersection).
   */
  GEOGen::EdgeFlags flags_;
};
}  // namespace matfp
