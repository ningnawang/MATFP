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
#include <geogram/basic/geometry.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/memory.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/mesh/mesh_sampling.h>

#include <algorithm>
#include <string>

#include "matfp/AABBWrapper.h"
#include "matfp/Logger.h"
#include "matfp/Parameters.h"

namespace matfp {

inline void mesh_find_sample_adjacent_faces(
    const GEO::Mesh& mesh,
    std::map<int, std::unordered_set<int>>& mesh_sample_adj_faces) {
  mesh_sample_adj_faces.clear();
  for (int f = 0; f < mesh.facets.nb(); f++) {
    for (int lv = 0; lv < 3; lv++) {
      GEO::index_t v = mesh.facets.vertex(f, lv);
      mesh_sample_adj_faces[v].insert(f);
    }
  }
}

// Debug
// only add downsample_percentage of the surface face centroids
template <GEO::index_t DIM>
inline bool
mesh_generate_adaptive_samples_on_surface_with_normals_given_downsample(
    const GEO::Mesh& mesh, std::vector<double>& p, std::vector<double>& normals,
    const double& downsample_percentage = 1.0,
    GEO::signed_index_t facets_begin_in = -1,
    GEO::signed_index_t facets_end_in = -1) {
  // not triangles, but polygons
  // geo_assert(mesh.facets.are_simplices());
  geo_assert(mesh.vertices.dimension() >= DIM);
  geo_assert(mesh.facets.nb() > 0);

  GEO::index_t facets_begin = 0;
  GEO::index_t facets_end = mesh.facets.nb();
  if (facets_begin_in != -1) {
    facets_begin = GEO::index_t(facets_begin_in);
  }
  if (facets_end_in != -1) {
    facets_end = GEO::index_t(facets_end_in);
  }

  typedef GEO::vecng<DIM, double> Point;
  GEO::signed_index_t first_t = -1;
  GEO::signed_index_t last_t = 0;

  for (GEO::index_t cur_t = facets_begin; cur_t < facets_end; cur_t++) {
    if (first_t == -1) {
      first_t = GEO::signed_index_t(cur_t);
    }
    last_t = std::max(last_t, GEO::signed_index_t(cur_t));
    GEO::index_t facet_nb_vertices = mesh.facets.nb_vertices(cur_t);
    GEO::index_t v1, v2, v3;  // split polygon to triangles, if not triangles
    for (GEO::index_t i = 1; i + 1 < facet_nb_vertices; i++) {
      // random number
      double r = ((double)std::rand() / (RAND_MAX));
      if (r > downsample_percentage) continue;

      // For debug: add only centroid seed per triangle
      GEO::vec3 cur_p = GEO::Geom::mesh_facet_center(mesh, cur_t);
      GEO::vec3 n = GEO::normalize(GEO::Geom::mesh_facet_normal(mesh, cur_t));
      for (GEO::coord_index_t coord = 0; coord < DIM; coord++) {
        p.push_back(cur_p[coord]);
        normals.push_back(n[coord]);
      }
    }  // for facet_nb_vertices done
  }    // for cur_t done

  if (mesh.facets.nb() > 1 && last_t == first_t) {
    GEO::Logger::warn("Sampler")
        << "Did put all the points in the same triangle" << std::endl;
    return false;
  }
  return true;
}

template <GEO::index_t DIM>
inline bool
mesh_generate_adaptive_samples_on_surface_with_normals_given_sizing_field(
    const GEO::Mesh& mesh, std::vector<double>& p, std::vector<double>& normals,
    GEO::Attribute<double>& densities, GEO::signed_index_t facets_begin_in = -1,
    GEO::signed_index_t facets_end_in = -1) {
  geo_assert(mesh.facets.are_simplices());
  geo_assert(mesh.vertices.dimension() >= DIM);
  geo_assert(mesh.facets.nb() > 0);

  GEO::index_t facets_begin = 0;
  GEO::index_t facets_end = mesh.facets.nb();
  if (facets_begin_in != -1) {
    facets_begin = GEO::index_t(facets_begin_in);
  }
  if (facets_end_in != -1) {
    facets_end = GEO::index_t(facets_end_in);
  }

  typedef GEO::vecng<DIM, double> Point;
  GEO::signed_index_t first_t = -1;
  GEO::signed_index_t last_t = 0;

  // For each triangle, we choose maximum weight/density maxWeight to
  // represent the whole triangle,
  // density rho = 1/(sizing field)^4 = 1/(r * LFS)^4
  // therefore
  // #estimated samples = Area(triangle) * sqrt(maxWeight) / (3 * sqrt(3) / 2.)
  double nb_seeds_per_facet = 0., fractpart, intpart;
  double cur_area = 0;
  for (GEO::index_t cur_t = facets_begin; cur_t < facets_end; cur_t++) {
    if (first_t == -1) {
      first_t = GEO::signed_index_t(cur_t);
    }
    last_t = std::max(last_t, GEO::signed_index_t(cur_t));

    double max_rho = DBL_MIN;
    // fetch minimum sizing field (max rho)
    for (GEO::coord_index_t coord = 0; coord < DIM; coord++) {
      max_rho = std::max(max_rho, densities[mesh.facets.vertex(cur_t, coord)]);
    }

    GEO::index_t v1 = mesh.facets.vertex(cur_t, 0);
    GEO::index_t v2 = mesh.facets.vertex(cur_t, 1);
    GEO::index_t v3 = mesh.facets.vertex(cur_t, 2);
    cur_area = GEO::Geom::triangle_area(
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v1)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v2)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v3)));
    nb_seeds_per_facet =
        cur_area * std::sqrt(std::fabs(max_rho)) / (3. * std::sqrt(3) / 2.);
    // nb_seeds_per_facet = cur_area * (std::fabs(max_rho)) / (3. * std::sqrt(3)
    // / 2.); add probablity for deciding floating part
    fractpart = std::modf(nb_seeds_per_facet, &intpart);
    double rand = ((double)std::rand() / (RAND_MAX));
    if (rand <= fractpart) {
      // logger().debug("nb_seeds_per_facet: {}, rand {} < fractpart {}, use
      // ceil",
      //     nb_seeds_per_facet, rand, fractpart
      // );
      nb_seeds_per_facet = std::ceil(nb_seeds_per_facet);
    } else {
      nb_seeds_per_facet = std::floor(nb_seeds_per_facet);
    }

    // add random seeds
    for (GEO::index_t add_cnt = 0; add_cnt < nb_seeds_per_facet; add_cnt++) {
      Point cur_p = GEO::Geom::random_point_in_triangle(
          *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v1)),
          *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v2)),
          *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v3)));
      GEO::vec3 n = GEO::normalize(GEO::Geom::mesh_facet_normal(mesh, cur_t));
      for (GEO::coord_index_t coord = 0; coord < DIM; coord++) {
        p.push_back(cur_p[coord]);
        normals.push_back(n[coord]);
      }
    }
  }  // for cur_t done

  if (mesh.facets.nb() > 1 && last_t == first_t) {
    GEO::Logger::warn("Sampler")
        << "Did put all the points in the same triangle" << std::endl;
    return false;
  }
  return true;
}

/**
 * \brief Generates a set of random samples over a surfacic mesh.
 * \param[in] mesh the mesh
 * \param[out] p pointer to an array of generated samples, of size
 *   \p nb_points times DIM. To be allocated by the caller.
 * \param[out] normals pointer to an array of facet normals corresponding to
 * generated samples, of size \p nb_points times DIM. To be allocated by the
 * caller. \param[in] nb_points number of points to generate \param[in] weight a
 * reference to a vertex attribute. If bound, it is taken into account.
 * \param[in] facets_begin_in if specified, first index of the facet
 *  sequence in which points should be generated. If left unspecified (-1),
 *  points are generated over all the facets of the mesh.
 * \param[in] facets_end_in if specified, one position past the last
 *  index of the facet sequence in which points should be generated.
 *  If left unspecified (-1), points are generated over all the facets
 *  of the mesh.
 * \tparam DIM dimension of the points, specified as a template argument
 *  for efficiency reasons
 * \return true if everything went OK, false otherwise. Whenever all the
 *  points land in the same facet, the function returns false to notify
 *  a potential numerical problem.
 */
template <GEO::index_t DIM>
inline bool mesh_generate_random_samples_on_surface_with_normals(
    const GEO::Mesh& mesh, double* p, double* normals, GEO::index_t nb_points,
    GEO::Attribute<double>& weight, GEO::signed_index_t facets_begin_in = -1,
    GEO::signed_index_t facets_end_in = -1) {
  geo_assert(mesh.facets.are_simplices());
  geo_assert(mesh.vertices.dimension() >= DIM);
  geo_assert(mesh.facets.nb() > 0);

  GEO::index_t facets_begin = 0;
  GEO::index_t facets_end = mesh.facets.nb();
  if (facets_begin_in != -1) {
    facets_begin = GEO::index_t(facets_begin_in);
  }
  if (facets_end_in != -1) {
    facets_end = GEO::index_t(facets_end_in);
  }

  typedef GEO::vecng<DIM, double> Point;

  // To ensure reproducibility accros successive
  // runs, reset the random number generator.
  GEO::Numeric::random_reset();

  GEO::vector<double> s(nb_points);
  for (GEO::index_t i = 0; i < nb_points; i++) {
    s[i] = GEO::Numeric::random_float64();
  }
  std::sort(s.begin(), s.end());

  double Atot = 0.0;
  for (GEO::index_t t = facets_begin; t < facets_end; ++t) {
    double At = GEO::mesh_facet_mass<DIM>(mesh, t, weight);
    Atot += At;
  }

  GEO::signed_index_t first_t = -1;
  GEO::signed_index_t last_t = 0;

  GEO::index_t cur_t = facets_begin;
  double cur_s = GEO::mesh_facet_mass<DIM>(mesh, facets_begin, weight) / Atot;
  for (GEO::index_t i = 0; i < nb_points; i++) {
    geo_debug_assert(i < s.size());
    while (s[i] > cur_s && cur_t < facets_end - 1) {
      cur_t++;
      geo_debug_assert(cur_t < facets_end);
      cur_s += GEO::mesh_facet_mass<DIM>(mesh, cur_t, weight) / Atot;
    }
    if (first_t == -1) {
      first_t = GEO::signed_index_t(cur_t);
    }
    last_t = std::max(last_t, GEO::signed_index_t(cur_t));

    // TODO: take weights into account
    //  with GEO::a new random_point_in_triangle_weighted()
    //  function.
    GEO::index_t v1 = mesh.facets.vertex(cur_t, 0);
    GEO::index_t v2 = mesh.facets.vertex(cur_t, 1);
    GEO::index_t v3 = mesh.facets.vertex(cur_t, 2);
    Point cur_p = GEO::Geom::random_point_in_triangle(
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v1)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v2)),
        *reinterpret_cast<const Point*>(mesh.vertices.point_ptr(v3)));
    GEO::vec3 n = GEO::normalize(GEO::Geom::mesh_facet_normal(mesh, cur_t));
    for (GEO::coord_index_t coord = 0; coord < DIM; coord++) {
      p[i * DIM + coord] = cur_p[coord];
      normals[i * DIM + coord] = n[coord];
    }
  }
  if (mesh.facets.nb() > 1 && last_t == first_t) {
    GEO::Logger::warn("Sampler")
        << "Did put all the points in the same triangle" << std::endl;
    return false;
  }
  return true;
}

}  // namespace matfp
