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
#include <matfp/Types/RegularTriangulation.h>
#include <matfp/geogram/generic_RPD_vertex.h>

namespace matfp {
// using namespace GEOGen;

/**
 * \brief Internal representation of polygons for GenericVoronoiDiagram.
 * \details Stores both geometrical and symbolic representations.
 * \note This is an internal implementation class used by
 *  GEO::RestrictedVoronoiDiagram. It is not meant to be
 *  used directly by client code.
 */
class PolygonCGAL {
 public:
  // /**
  //  * \brief Internal representation of vertices.
  //  */
  // typedef matfp::Vertex Vertex;

  /**
   * \brief Gets the number of vertices.
   */
  index_t nb_vertices() const { return index_t(vertex_.size()); }

  /**
   * \brief Gets a vertex by index.
   * \param[in] i index of the Vertex in this Polygon
   * \return a const reference to the Vertex at index \p i
   * \pre \p i < nb_vertices()
   */
  const matfp::Vertex& vertex(index_t i) const {
    geo_debug_assert(i < nb_vertices());
    return vertex_[i];
  }

  /**
   * \brief Gets a vertex by index.
   * \param[in] i index of the Vertex in this Polygon
   * \return a reference to the Vertex at index \p i
   * \pre \p i < nb_vertices()
   */
  matfp::Vertex& vertex(index_t i) {
    geo_debug_assert(i < nb_vertices());
    return vertex_[i];
  }

  /**
   * \brief Gets the index of the successor of a Vertex.
   * \param[in] i index of the Vertex in this Polygon
   * \return the index of the successor of Vertex \p i
   * \pre \p i < nb_vertices()
   */
  index_t next_vertex(index_t i) const {
    geo_debug_assert(i < nb_vertices());
    return (i == nb_vertices() - 1) ? 0 : (i + 1);
  }

  /**
   * \brief Gets the index of the predecessor of a Vertex.
   * \param[in] i index of the Vertex in this Polygon
   * \return the index of the predecessor of Vertex \p i
   * \pre \p ii < nb_vertices()
   */
  index_t prev_vertex(index_t i) const {
    geo_debug_assert(i < nb_vertices());
    return (i == 0) ? (nb_vertices() - 1) : (i - 1);
  }

  /**
   * \brief Adds a Vertex to this Polygon.
   * \param[in] v the vertex to be added. It is copied.
   * \return the address of the stored vertex.
   */
  matfp::Vertex* add_vertex(const matfp::Vertex& v) {
    vertex_.push_back(v);
    return &*(vertex_.rbegin());
  }

  /**
   * \brief Clears this Polygon.
   */
  void clear() { vertex_.resize(0); }

  /**
   * \brief Resizes this Polygon.
   * \param[in] sz new size
   */
  void resize(index_t sz) { vertex_.resize(sz); }

  /**
   * \brief Assigns a mesh facet to this Polygon.
   * \details The facet from the initial mesh is converted into
   *  the internal geometric/symbolic representation.
   * \param[in] mesh the mesh from which the facet is copied
   * \param[in] f the index of the facet in \p mesh
   * \param[in] symbolic if true, symbolic information is copied
   * \param[in] vertex_weight a reference to a vertex attribute
   *  that stores weights. If not bound, then 1.0 is used for
   *  the weights.
   */
  void initialize_from_mesh_facet(const GEO::Mesh* mesh, index_t f,
                                  bool symbolic,
                                  const GEO::Attribute<double>& vertex_weight);

  /**
   * \brief Clips a polygon with a plane.
   * \details Computes the intersection between this Polygon
   * and the half-space determined by the positive side
   * of the bisector of segment [i,j] (on the same side as vertex i).
   *
   * \param[out] target where to store the intersection
   * \param[out] target_intersections
   *    where to allocate the generated vertices
   * \param[in] mesh the input mesh, used by the symbolic information
   * \param[in] delaunay the Delaunay triangulation
   * \param[in] i index of one extremity of bisector in \p delaunay
   * \param[in] j index of the other extremity of the bisector
   *    in \p delaunay
   * \param[in] exact if true, exact predicates are used.
   *   Implies symbolic.
   * \param[in] symbolic if true, symbolic representation
   *   of vertices is computed
   */
  template <index_t DIM>
  void clip_by_plane(PolygonCGAL& target, PointAllocator& target_intersections,
                     const GEO::Mesh* mesh, const RegularTriangulationNN* rt,
                     Vertex_handle_rt& i, Vertex_handle_rt& j, bool exact,
                     bool symbolic) {
    // if(exact) {
    clip_by_plane_exact<3>(target, target_intersections, mesh, rt, i, j);
    // } else {
    //     clip_by_plane_fast<DIM>(
    //         target, target_intersections, rt, i, j, symbolic
    //     );
    // }
  }

  /**
   * \brief Overwrites this Polygon with the contents of another
   *  polygon.
   * \param[in] rhs a const reference to the polygon to be copied.
   */
  void copy(const PolygonCGAL& rhs) { vertex_ = rhs.vertex_; }

  /**
   * \brief Swaps the contents of this Polygon and another polygon.
   * \param[in,out] rhs a reference to the Polygon to be swapped with
   *  this one.
   */
  void swap(PolygonCGAL& rhs) { vertex_.swap(rhs.vertex_); }

 protected:
  /**
   * \brief Clips a Polygon with a plane (fast inexact version).
   * \details Computes the intersection between this Polygon
   * and the half-space determined by the positive side
   * of the bisector of segment [i,j] (the side of i).
   * This version uses a "fused" predicates-constructions
   * strategy (and reuses the computations from the predicates
   * to accelerate the constructions).
   *
   * \param[out] target where to store the intersection
   * \param[out] target_intersections
   *   where to allocate the generated vertices
   * \param[in] delaunay the Delaunay triangulation
   * \param[in] i index of one extremity of bisector in \p delaunay
   * \param[in] j index of the other extremity
   *   of the bisector in \p delaunay
   * \param[in] symbolic if true, symbolic representation
   *  of vertices is computed
   *
   * \internal
   * \note Profiling revealed that this routine is where
   * the system spends the largest amount of time
   * (no big surprise...).
   */
  template <index_t DIM>
  void clip_by_plane_fast(PolygonCGAL& target,
                          PointAllocator& target_intersections,
                          const RegularTriangulationNN* rt, Vertex_handle_rt i,
                          Vertex_handle_rt j, bool symbolic) const {
    target.clear();
    if (nb_vertices() == 0) {
      return;
    }

    // TODO: remember to delete if using this function
    const double* geo_restrict pi = rt->get_double_data(i->point());
    geo_assume_aligned(pi, geo_dim_alignment(DIM));
    const double* geo_restrict pj = rt->get_double_data(j->point());
    geo_assume_aligned(pj, geo_dim_alignment(DIM));

    // Compute d = n . m, where n is the
    // normal vector of the bisector [pi,pj]
    // and m the middle point of the bisector.
    geo_decl_aligned(double d);
    d = 0;
    for (coord_index_t c = 0; c < DIM; ++c) {
      d += (pi[c] + pj[c]) * (pi[c] - pj[c]);
    }

    // The predecessor of the first vertex is the last vertex
    index_t prev_k = nb_vertices() - 1;
    const Vertex* prev_vk = &(vertex(prev_k));
    const double* geo_restrict prev_pk = prev_vk->point();
    geo_assume_aligned(prev_pk, geo_dim_alignment(DIM));

    // We compute:
    //    prev_l = prev_vk . n
    geo_decl_aligned(double prev_l);
    prev_l = 0.0;
    for (coord_index_t c = 0; c < DIM; ++c) {
      prev_l += prev_pk[c] * (pi[c] - pj[c]);
    }

    // We compute:
    //    side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
    GEO::Sign prev_status = GEO::geo_sgn(2.0 * prev_l - d);

    for (index_t k = 0; k < nb_vertices(); k++) {
      const Vertex* vk = &(vertex(k));
      const double* pk = vk->point();

      // We compute: l = vk . n
      geo_decl_aligned(double l);
      l = 0.0;
      for (coord_index_t c = 0; c < DIM; ++c) {
        l += pk[c] * (pi[c] - pj[c]);
      }

      // We compute:
      //   side1(pi,pj,q) = sign(2*q.n - n.m) = sign(2*l - d)
      GEO::Sign status = GEO::geo_sgn(2.0 * l - d);

      // If status of edge extremities differ,
      // then there is an intersection.
      if (status != prev_status && (prev_status != 0)) {
        Vertex I;
        double* Ipoint = target_intersections.new_item();
        I.set_point(Ipoint);
        if (symbolic) {
          if (!I.sym().intersect_symbolic(prev_vk->sym(), vk->sym(),
                                          j->info().tag)) {
            // We encountered a problem. As a workaround,
            // we copy prev_vk into the result.
            I = *prev_vk;
          }
        }

        // Compute lambda1 and lambda2, the
        // barycentric coordinates of the intersection I
        // in the segment [prev_vk vk]
        // Note that d and l (used for the predicates)
        // are reused here.
        double denom = 2.0 * (prev_l - l);
        double lambda1, lambda2;

        // Shit happens ! [Forrest Gump]
        if (::fabs(denom) < 1e-20) {
          lambda1 = 0.5;
          lambda2 = 0.5;
        } else {
          lambda1 = (d - 2.0 * l) / denom;
          // Note: lambda2 is also given
          // by (2.0*l2-d)/denom
          // (but 1.0 - lambda1 is a bit
          //  faster to compute...)
          lambda2 = 1.0 - lambda1;
        }
        // Compute intersection I by weighting
        // the edge extremities with the barycentric
        // coordinates lambda1 and lambda2
        for (coord_index_t c = 0; c < DIM; ++c) {
          Ipoint[c] = lambda1 * prev_pk[c] + lambda2 * pk[c];
        }
        I.set_weight(lambda1 * prev_vk->weight() + lambda2 * vk->weight());
        if (status > 0) {
          I.copy_edge_from(*prev_vk);
          I.set_adjacent_seed(signed_index_t(j->info().tag));
        } else {
          I.set_flag(GEOGen::INTERSECT);
          I.set_adjacent_seed(vk->adjacent_seed());
        }
        target.add_vertex(I);
      }
      if (status > 0) {
        target.add_vertex(*vk);
      }
      prev_vk = vk;
      prev_pk = pk;
      prev_status = status;
      prev_k = k;
      prev_l = l;
    }
  }

  /**
   * \brief Clips a Polygon with a plane (exact version).
   * \details Computes the intersection between this Polygon
   * and the half-space determined by the positive side
   * of the bisector of segment [i,j] (the side of i).
   * This version uses symbolically perturbed exact predicates.
   *
   * \param[out] target where to store the intersection
   * \param[out] target_intersections
   *  where to allocate the generated vertices
   * \param[in] mesh the input mesh (used by exact predicates)
   * \param[in] delaunay the Delaunay triangulation
   * \param[in] i index of one extremity of bisector in \p delaunay
   * \param[in] j index of the other extremity of
   *  the bisector in \p delaunay
   */
  template <index_t DIM>
  void clip_by_plane_exact(PolygonCGAL& target,
                           PointAllocator& target_intersections,
                           const GEO::Mesh* mesh,
                           const RegularTriangulationNN* rt,
                           Vertex_handle_rt& i, Vertex_handle_rt& j) {
    // logger().debug("In clip_by_plane_exact-----");
    // logger().debug("nb_vertices {}", nb_vertices());

    // std::cout << "-- In clip_by_plane_exact --" << std::endl;
    // std::cout << "nb_vertices: " << nb_vertices() << std::endl;

    target.clear();
    if (nb_vertices() == 0) {
      return;
    }

    // logger().debug("i pos ({},{},{}), j pos ({},{},{})",
    //     ip[0], ip[1], ip[2], jp[0], jp[1], jp[2]
    // );

    // std::vector<double> pii{ip[0], ip[1], ip[2]};
    // std::vector<double> pjj{jp[0], jp[1], jp[2]};

    // const double* pi = pii.data();
    // const double* pj = pjj.data();

    const double* pi = rt->get_double_data(i->point());
    const double* pj = rt->get_double_data(j->point());

    const double wi = rt->get_weight(i->point());
    const double wj = rt->get_weight(j->point());

    // if (i->info().tag == 457) {
    //     logger().debug("pi pos {}: ({},{},{},{}) , pj pos {}: ({},{},{},{})",
    //         i->info().tag, pi[0], pi[1], pi[2], pi[3],
    //         j->info().tag, pj[0], pj[1], pj[2], pj[3]
    //     );
    // }

    // std::cout << "--- processisng pi: " << i->info().tag << " and pj: " <<
    // j->info().tag << std::endl;

    // The predecessor of the first vertex is the last vertex
    index_t prev_k = nb_vertices() - 1;
    const Vertex* prev_vk = &(vertex(prev_k));
    GEO::Sign prev_status = side_exact(mesh, rt, *prev_vk, pi, wi, pj, wj, DIM);

    const double* prev_vkp = prev_vk->point();

    // if (i->info().tag == 457) {
    //     print_sign(prev_status);
    //     logger().debug("prev_vk: ({},{},{})",
    //         prev_vkp[0], prev_vkp[1], prev_vkp[2]
    //     );
    // }

    // print_sign(prev_status);
    // std::cout << "prev_vkp: (" << prev_vkp[0] << ", " << prev_vkp[1] << ", "
    // << prev_vkp[2] << ")" << std::endl;

    for (index_t k = 0; k < nb_vertices(); ++k) {
      const matfp::Vertex* vk = &(vertex(k));
      const double* vkp = vk->point();

      // if (i->info().tag == 457) {
      //     logger().debug("vkp: ({},{},{})",
      //         vkp[0], vkp[1], vkp[2]
      //     );
      // }

      GEO::Sign status = side_exact(mesh, rt, *vk, pi, wi, pj, wj, DIM);

      // if (i->info().tag == 457) {
      // print_sign(status);
      // std::cout << "vkp: (" << vkp[0] << ", " << vkp[1] << ", " << vkp[2] <<
      // ")" << std::endl;
      // }

      // If status of edge extremities differ,
      // there is an intersection.
      if (status != prev_status && (prev_status != 0)) {
        matfp::Vertex I;
        if (!I.sym().intersect_symbolic(prev_vk->sym(), vk->sym(),
                                        j->info().tag)) {
          // We encountered a problem. As a workaround,
          // we copy prev_vk into the result.
          I = *prev_vk;
          // geo_assert_not_reached ;
          // not supposed to happen in exact mode
          // logger().error("++++++++++++++++ we have a issue here with
          // intersect_symbolic!!!");
        }
        I.intersect_geom<DIM>(target_intersections, *prev_vk, *vk, pi, pj);

        // I.intersect_geom_new<DIM>(
        //     target_intersections, *prev_vk, *vk, pi, pj
        // );

        // logger().debug("created new point II ({},{},{})", I.point()[0],
        // I.point()[1], I.point()[2]);
        if (status > 0) {
          I.copy_edge_from(*prev_vk);
          I.set_adjacent_seed(signed_index_t(j->info().tag));
          // std::cout << "status > 0, new point I has adjacent_seed: " <<
          // I.adjacent_seed() << std::endl;
        } else {
          I.set_flag(GEOGen::INTERSECT);
          I.set_adjacent_seed(vk->adjacent_seed());
          // std::cout << "status <= 0, new point I has adjacent_seed: " <<
          // I.adjacent_seed() << std::endl;
        }
        target.add_vertex(I);

        // if (i->info().tag == 457) {
        //     logger().debug("created new point I ({},{},{})", I.point()[0],
        //     I.point()[1], I.point()[2]); logger().debug("status > 0 ? {},
        //     I.adjacent_seed(): {}", status > 0, I.adjacent_seed());
        //     logger().debug("I.sym().nb_bisectors(): {}",
        //     I.sym().nb_bisectors()); for (int i = 0; i <
        //     I.sym().nb_bisectors(); i++) {
        //         logger().debug("I has bisector {}: {}", i,
        //         I.sym().bisector(i));
        //     }
        // }

        // std::cout << "created new point I (" << I.point()[0] << ", " <<
        // I.point()[1] << ", " << I.point()[2] << ")" << std::endl;
      }
      if (status > 0) {
        target.add_vertex(*vk);
      }
      prev_vk = vk;
      prev_status = status;
      prev_k = k;
    }
    delete[] pi;
    delete[] pj;
    pi = nullptr;
    pj = nullptr;

    // if (i->info().tag == 457) {
    //     logger().debug("target nb_vertices {}", target.nb_vertices());
    // }

    // std::cout << "target nb_vertices: " << target.nb_vertices() << std::endl;
  }

  inline void print_sign(GEO::Sign& s) {
    switch (s) {
      case GEO::Sign::ZERO:
        // logger().debug("sign is 0");
        std::cout << "sign is 0" << std::endl;
        break;
      case GEO::Sign::POSITIVE:
        // logger().debug("sign is 1");
        std::cout << "sign is 1" << std::endl;
        break;
      case GEO::Sign::NEGATIVE:
        // logger().debug("sign is -1");
        std::cout << "sign is -1" << std::endl;
        break;
      default:
        break;
    }
  }

  /**
   * \brief Returns the position of a point
   * relative to a bisector (exact version).
   * \details Position of q relative to the bisector Pi(i,j).
   *  The symbolic representation of q is used. Symbolic
   *  perturbation is applied to degenerate configurations,
   *  therefore ZERO is never returned.
   * \param[in] mesh the input mesh
   * \param[in] delaunay the Delaunay triangulation
   * \param[in] q query point
   * \param[in] pi one extremity of the bisector
   * \param[in] pj the other extremity of the bisector
   * \param[in] dim dimension of the points
   * \return POSITIVE if q is on pi's side, NEGATIVE otherwise
   *  (ZERO is never encountered thanks to globally coherent
   *  symbolic perturbations).
   */
  static GEO::Sign side_exact(const GEO::Mesh* mesh,
                              const RegularTriangulationNN* rt,
                              const matfp::Vertex& q, const double* pi,
                              double wi, const double* pj, double wj,
                              GEO::coord_index_t dim);

 private:
  GEO::vector<matfp::Vertex> vertex_;
};

}  // namespace matfp