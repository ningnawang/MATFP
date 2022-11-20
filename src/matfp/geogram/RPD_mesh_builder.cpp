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

#include <matfp/Logger.h>
#include <matfp/geogram/RPD_mesh_builder.h>

namespace matfp {

RPDVertexMap::RPDVertexMap() : nb_vertices_(0) {}

GEO::index_t RPDVertexMap::find_or_create_vertex(
    GEO::index_t center_vertex_id, const matfp::SymbolicVertex& sym) {
  // if (center_vertex_id == 4 || center_vertex_id == 5) {
  //     logger().debug("center_vertex_id: {}, sym.nb_bisectors: {}",
  //         center_vertex_id, sym.nb_bisectors()
  //     );
  // }

  switch (sym.nb_bisectors()) {
    case 3: {
      GEO::index_t ib1 = sym.bisector(0);
      GEO::index_t ib2 = sym.bisector(1);
      GEO::index_t ib3 = sym.bisector(2);
      quadindex K(center_vertex_id + 1, ib1 + 1, ib2 + 1, ib3 + 1);
      auto it = ppp_to_id_.find(K);
      if (it != ppp_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();

        if (result == 15477 || result == 15555) {
          logger().debug(
              "vertex with id {}: sym.nb_bisectors {}, ib1 {}, ib2 {}, ib3 {}",
              result, sym.nb_bisectors(), ib1, ib2, ib3);
        }

        ppp_to_id_[K] = result;
        return result;
      }
    }
    case 2: {
      GEO::index_t f = sym.boundary_facet(0);
      GEO::index_t ib1 = sym.bisector(0);
      GEO::index_t ib2 = sym.bisector(1);
      signed_quadindex K(GEO::signed_index_t(center_vertex_id) + 1,
                         -GEO::signed_index_t(f) - 1,
                         GEO::signed_index_t(ib1) + 1,
                         GEO::signed_index_t(ib2) + 1);

      // if (center_vertex_id == 4 || center_vertex_id == 5) {
      //     logger().debug("center_vertex_id: {},
      //     GEO::signed_index_t(center_vertex_id)+1: {}",
      //         center_vertex_id, GEO::signed_index_t(center_vertex_id)+1
      //     );
      //     logger().debug(" GEO::signed_index_t(ib1)+1: {},
      //     GEO::signed_index_t(ib2)+1: {}",
      //         GEO::signed_index_t(ib1) + 1,  GEO::signed_index_t(ib2) + 1
      //     );
      //     logger().debug("signed_quadindex K ({},{},{},{})",
      //         K.indices[0], K.indices[1], K.indices[2], K.indices[3]
      //     );
      // }

      auto it = ppm_to_id_.find(K);
      if (it != ppm_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();

        if (result == 15477 || result == 15555) {
          logger().debug(
              "vertex with id {}: sym.nb_bisectors {}, ib1 {}, ib2 {}", result,
              sym.nb_bisectors(), ib1, ib2);
        }

        ppm_to_id_[K] = result;
        return result;
      }
    }
    case 1: {
      GEO::index_t bv1, bv2;
      sym.get_boundary_edge(bv1, bv2);
      GEO::index_t ib = sym.bisector(0);
      signed_quadindex K(GEO::signed_index_t(center_vertex_id) + 1,
                         -GEO::signed_index_t(bv1) - 1,
                         -GEO::signed_index_t(bv2) - 1,
                         GEO::signed_index_t(ib) + 1);
      auto it = pmm_to_id_.find(K);
      if (it != pmm_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();

        if (result == 15477 || result == 15555) {
          logger().debug("vertex with id {}: sym.nb_bisectors {}, ib1 {}",
                         result, sym.nb_bisectors(), ib);
        }

        pmm_to_id_[K] = result;
        return result;
      }
    }
    case 0: {
      GEO::index_t bv = sym.get_boundary_vertex();
      if (bv >= bv_to_id_.size()) {
        bv_to_id_.resize(bv + 1, -1);
      }
      if (bv_to_id_[bv] == -1) {
        bv_to_id_[bv] = GEO::signed_index_t(new_vertex());
      }

      if (bv_to_id_[bv] == 15477 || bv_to_id_[bv] == 15555) {
        logger().debug("vertex with id {}: sym.nb_bisectors {}", bv_to_id_[bv],
                       sym.nb_bisectors());
      }

      return GEO::index_t(bv_to_id_[bv]);
    }
    default:
      geo_assert_not_reached;
  }
}

GEO::index_t RPDVertexMap::find_or_create_vertex(
    GEO::index_t center_vertex_id, const matfp::SymbolicVertex& sym,
    std::set<GEO::index_t>& v_adj) {
  v_adj.clear();  // seems okay?
  switch (sym.nb_bisectors()) {
    case 3: {
      GEO::index_t ib1 = sym.bisector(0);
      GEO::index_t ib2 = sym.bisector(1);
      GEO::index_t ib3 = sym.bisector(2);
      v_adj.insert(ib1);
      v_adj.insert(ib2);
      v_adj.insert(ib3);

      quadindex K(center_vertex_id + 1, ib1 + 1, ib2 + 1, ib3 + 1);
      auto it = ppp_to_id_.find(K);
      if (it != ppp_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();
        ppp_to_id_[K] = result;
        return result;
      }
    }
    case 2: {
      GEO::index_t f = sym.boundary_facet(0);
      GEO::index_t ib1 = sym.bisector(0);
      GEO::index_t ib2 = sym.bisector(1);
      v_adj.insert(ib1);
      v_adj.insert(ib2);
      signed_quadindex K(GEO::signed_index_t(center_vertex_id) + 1,
                         -GEO::signed_index_t(f) - 1,
                         GEO::signed_index_t(ib1) + 1,
                         GEO::signed_index_t(ib2) + 1);

      auto it = ppm_to_id_.find(K);
      if (it != ppm_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();
        ppm_to_id_[K] = result;
        return result;
      }
    }
    case 1: {
      GEO::index_t bv1, bv2;
      sym.get_boundary_edge(bv1, bv2);
      GEO::index_t ib = sym.bisector(0);
      v_adj.insert(ib);

      signed_quadindex K(GEO::signed_index_t(center_vertex_id) + 1,
                         -GEO::signed_index_t(bv1) - 1,
                         -GEO::signed_index_t(bv2) - 1,
                         GEO::signed_index_t(ib) + 1);
      auto it = pmm_to_id_.find(K);
      if (it != pmm_to_id_.end()) {
        return it->second;
      } else {
        GEO::index_t result = new_vertex();
        pmm_to_id_[K] = result;
        return result;
      }
    }
    case 0: {
      GEO::index_t bv = sym.get_boundary_vertex();
      if (bv >= bv_to_id_.size()) {
        bv_to_id_.resize(bv + 1, -1);
      }
      if (bv_to_id_[bv] == -1) {
        bv_to_id_[bv] = GEO::signed_index_t(new_vertex());
      }
      return GEO::index_t(bv_to_id_[bv]);
    }
    default:
      geo_assert_not_reached;
  }
}
}  // namespace matfp
