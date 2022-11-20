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

#include <geogram/basic/common.h>
#include <geogram/basic/memory.h>
#include <matfp/Types/RegularTriangulation.h>
#include <matfp/geogram/generic_RPD_cell.h>
#include <matfp/geogram/generic_RPD_polygon.h>
#include <matfp/geogram/generic_RPD_vertex.h>

#include <stack>
#include <vector>

namespace matfp {
using namespace GEO;
/**
 * \brief A (facet,seed) pair.
 * \details Used by matfp::RestrictedPowerDiagram
 * for propagating over the facet graph and the
 * Delaunay 1-skeleton.
 */
struct FacetSeedHandle {
  /**
   * \brief Creates a new FacetSeedHandle
   * \param[in] f_in index of the facet
   * \param[in] seed_in Vertex_handle of the seed
   */
  FacetSeedHandle(GEO::index_t f_in, Vertex_handle_rt seed_in)
      : f(f_in), seed(seed_in) {}

  /**
   * \brief Creates a new uninitialized FacetSeedHandle
   * \details F and seed contain random values.
   */
  FacetSeedHandle() {}

  /**
   * \brief Compares two facet seeds using lexicographic order.
   * \details Makes it possible to use FacetSeedHandle as keys for
   * std::set and std::map.
   */
  bool operator<(const FacetSeedHandle& rhs) const {
    if (f < rhs.f) {
      return true;
    }
    if (f > rhs.f) {
      return false;
    }
    return seed->info().tag < rhs.seed->info().tag;
  }

  GEO::index_t f;
  Vertex_handle_rt seed;
};

/**
 * \brief A (tetrahedron,seed) pair.
 * \details Used by GEOGen::RestrictedVoronoiDiagram
 * for propagating over the tetrahedra graph and the
 * Delaunay 1-skeleton.
 */
typedef FacetSeedHandle TetSeedHandle;

/************************************************************************/

/**
 * \brief A stack of FacetSeed.
 * \details Used by GEOGen::RestrictedVoronoiDiagram.
 */
typedef std::stack<FacetSeedHandle> FacetSeedHandleStack;

/**
 * \brief A stack of TetSeed.
 * \details Used by GEOGen::RestrictedVoronoiDiagram.
 */
typedef std::stack<TetSeedHandle> TetSeedHandleStack;
/**
 * \brief A stack of seed indices (index_t).
 * \details Used by GEOGen::RestrictedVoronoiDiagram.
 */
typedef std::stack<Vertex_handle_rt> SeedHandleStack;

}  // namespace matfp