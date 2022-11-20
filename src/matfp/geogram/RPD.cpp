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
#include "matfp/geogram/RPD.h"

#include <geogram/basic/algorithm.h>
#include <geogram/basic/argused.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/geometry_nd.h>
#include <geogram/basic/process.h>
#include <geogram/bibliography/bibliography.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/mesh/mesh_partition.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_sampling.h>
#include <geogram/voronoi/integration_simplex.h>

#include "matfp/Common.h"
#include "matfp/Logger.h"
#include "matfp/geogram/RPD_callback.h"
#include "matfp/geogram/RPD_mesh_builder.h"
#include "matfp/geogram/generic_RPD.h"

namespace {
using namespace GEO;
using namespace matfp;

/**
 * \brief Generic implementation of RestrictedPowerDiagram.
 * \tparam DIM dimension
 */
class RPD_3d_Impl : public matfp::RestrictedPowerDiagram {
  /** \brief This class type */
  typedef RPD_3d_Impl thisclass;

  /** \brief The base class of this class */
  typedef RestrictedPowerDiagram baseclass;

 public:
  /** \brief Implementation based on the generic version. */
  typedef matfp::GenRestrictedPowerDiagram GenRestrictedPowerDiagram;

  /** \brief Representation of points. */
  typedef vecng<3, double> Point;

  /** \brief Representation of vectors. */
  typedef vecng<3, double> Vector;

  /** \brief Represents a point and its symbolic information. */
  typedef typename GenRestrictedPowerDiagram::Vertex Vertex;

  /**
   * \brief Specifies the computation done by the threads.
   */
  enum ThreadMode {
    MT_NONE,      /**< uninitialized                          */
    MT_LLOYD,     /**< Lloyd iteration                        */
    MT_NEWTON,    /**< Newton optimization                    */
    MT_INT_SMPLX, /**< Newton with integration simplex        */
    MT_POLYG,     /**< Polygon callback                       */
    MT_POLYH,     /**< Polyhedron callback                    */
    MT_RPD_S_MESH /**< Build surfacic RPD Mesh                */
  };

  /**
   * \brief Creates a RPD_3d_Impl.
   *
   * \details The dimension is determined by \p mesh->dimension().
   * \param[in] delaunay the Delaunay triangulation
   * \param[in] mesh the input mesh
   * \param[in] R3_embedding gives for each vertex
   *  its mapping in 3D space.
   * \param[in] R3_embedding_stride gives the stride between
   *  two consecutive vertices in R3_embedding
   */
  RPD_3d_Impl(RegularTriangulationNN* rt, Mesh* mesh,
              const double* R3_embedding, index_t R3_embedding_stride)
      : RestrictedPowerDiagram(rt, mesh, R3_embedding, R3_embedding_stride),
        RPD_(rt, mesh) {
    use_exact_projection_ = false;
    is_slave_ = false;
    master_ = nullptr;
    parts_ = nullptr;
    nb_parts_ = 0;
    funcval_ = 0.0;
    simplex_func_ = nullptr;
    polygon_callback_ = nullptr;
    polyhedron_callback_ = nullptr;
    build_rpd_ = nullptr;
    arg_vectors_ = nullptr;
    arg_scalars_ = nullptr;
    thread_mode_ = MT_NONE;
    nb_triangles_ = 0;
  }

  /**
   * \brief Constructor for parts, used in multithreading mode.
   */
  RPD_3d_Impl()
      : RestrictedPowerDiagram(nullptr, nullptr, nullptr, 0),
        RPD_(nullptr, nullptr) {
    use_exact_projection_ = false;
    is_slave_ = true;
    master_ = nullptr;
    mesh_ = nullptr;
    parts_ = nullptr;
    nb_parts_ = 0;
    facets_begin_ = -1;
    facets_end_ = -1;
    funcval_ = 0.0;
    simplex_func_ = nullptr;
    polygon_callback_ = nullptr;
    polyhedron_callback_ = nullptr;
    build_rpd_ = nullptr;
    arg_vectors_ = nullptr;
    arg_scalars_ = nullptr;
    thread_mode_ = MT_NONE;
    nb_triangles_ = 0;
  }

  /********************************************************************/

  /**
   * \brief Adapter class used internally to implement for_each_polygon()
   * \details Gets the current triangle from the RPD and passes it back
   *  to the callback. It is needed because GenericRPD::for_each_polygon()
   *  does not pass the current triangle.
   */
  // TODO: pass it through all the callbacks, because it is ridiculous:
  // we pass it through the first levels, then throw it, then retrieve it
  // (see GenRPD)
  class PolygonCallbackAction {
   public:
    /**
     * \brief PolygonCallbackAction constructor
     * \param[in] RPD a pointer to the restricted Power diagram
     * \param[in] callback a pointer to the PolygonCallback
     */
    PolygonCallbackAction(GenRestrictedPowerDiagram& RPD,
                          matfp::RPDPolygonCGALCallback& callback)
        : RPD_(RPD), callback_(callback) {}

    /**
     * \brief Callback called for each polygon.
     * \details Routes the callback to the wrapped user action class.
     * \param[in] v index of current Delaunay seed
     * \param[in] P intersection between current mesh facet
     *  and the Voronoi cell of \p v
     */
    void operator()(GEO::index_t v,
                    const GenRestrictedPowerDiagram::Polygon& P) const {
      callback_(v, RPD_.current_facet(), P);
    }

   protected:
    GenRestrictedPowerDiagram& RPD_;
    matfp::RPDPolygonCGALCallback& callback_;
  };

  /********************************************************************/

  /**
   * Helper function for computing with polygon callback in parallel
   **/
  virtual void compute_with_polygon_callback(
      matfp::RPDPolygonCGALCallback& polygon_callback) {
    logger().debug("calling compute_with_polygon_callback ...");
    create_threads();
    if (nb_parts() == 0) {
      logger().debug("in compute_with_polygon_callback nb_parts() == 0");
      PolygonCallbackAction action(RPD_, polygon_callback);
      RPD_.for_each_polygon(action);
    } else {
      logger().debug("in compute_with_polygon_callback nb_parts() != 0");

      for (index_t t = 0; t < nb_parts(); t++) {
        part(t).RPD_.set_symbolic(RPD_.symbolic());
        part(t).RPD_.set_connected_components_priority(
            RPD_.connected_components_priority());
      }
      spinlocks_.resize(rt_->get_nb_vertices());
      thread_mode_ = MT_POLYG;
      polygon_callback_ = &polygon_callback;
      polygon_callback_->set_spinlocks(&spinlocks_);
      // Note: callback begin()/end() is called in for_each_polygon()
      parallel_for(0, nb_parts(), [this](index_t i) { run_thread(i); });
      polygon_callback_->set_spinlocks(nullptr);
    }
  }

  virtual void compute_with_polyhedron_callback(
      matfp::RPDPolyhedronCGALCallback& polyhedron_callback) {
    create_threads();
    if (nb_parts() == 0) {
      RPD_.for_each_polyhedron(polyhedron_callback);
    } else {
      for (index_t t = 0; t < nb_parts(); t++) {
        part(t).RPD_.set_symbolic(RPD_.symbolic());
        part(t).RPD_.set_connected_components_priority(
            RPD_.connected_components_priority());
      }
      spinlocks_.resize(rt_->get_nb_vertices());
      thread_mode_ = MT_POLYH;
      polyhedron_callback_ = &polyhedron_callback;
      polyhedron_callback_->set_spinlocks(&spinlocks_);
      // Note: callback begin()/end() is
      // called in for_each_polyhedron()
      parallel_for(0, nb_parts(), [this](index_t i) { run_thread(i); });
      polyhedron_callback_->set_spinlocks(nullptr);
    }
  }

  /********************************************************************/
  /**
   * \brief Implementation class for explicitly constructing
   *    a surfacic mesh that corresponds to the surfacic
   *    restricted Voronoi diagram.
   * \details To be used as a template argument
   *    to RPD::for_each_polygon(). The current Vornoi cell is
   *    reported in facet region.
   * \tparam BUILDER a class that implements iterative mesh building,
   *    e.g., MeshBuilder.
   */
  template <class BUILDER>
  class BuildRPD {
   public:
    /**
     * \brief Constructs a new BuildRPD.
     * \param[in] RPD_in the restricted Voronoi diagram
     * \param[in] builder the lesh builder
     */
    BuildRPD(const GenRestrictedPowerDiagram& RPD_in, BUILDER& builder)
        : RPD(RPD_in), builder_(builder), current_facet_(-1) {
      // std::cout << "calling builder_.begin_surface() ...\n";
      builder_.begin_surface();
      global_lock_ = GEOGRAM_SPINLOCK_INIT;
      // std::cout << "done calling builder_.begin_surface() ...\n";
    }

    /**
     * \brief The destructor
     * \details Terminates the current facet
     *    and the current surface.
     */
    ~BuildRPD() {
      // std::cout << "calling ~BuildRPD() ..." << std::endl;
      if (current_facet_ != -1) {
        builder_.end_reference_facet();
      }
      builder_.end_surface();
      // std::cout << "finish ~BuildRPD()" << std::endl;
    }

    /**
     * \brief The callback called for each restricted Voronoi cell.
     * \param[in] v index of current center vertex
     * \param[in] P current restricted Voronoi cell
     */
    void operator()(GEO::index_t v,
                    const typename GenRestrictedPowerDiagram::Polygon& P) {
      GEO::Process::acquire_spinlock(global_lock_);
      // std::cout << "In BuildRPD operator ..." << std::endl;
      // std::cout << "processing seed: " << (int) v << std::endl;
      index_t f = RPD.current_facet();
      if (signed_index_t(f) != current_facet_) {
        if (current_facet_ != -1) {
          builder_.end_reference_facet();
        }
        current_facet_ = signed_index_t(f);
        builder_.begin_reference_facet(f);
      }

      // logger().debug("curren center vertex {}, nb_vertices {}", v,
      // P.nb_vertices()); std::cout << "builder start to build facet for seed "
      // << v << ", nb_vertices: " << P.nb_vertices() << std::endl;

      builder_.begin_facet(v);
      for (index_t i = 0; i < P.nb_vertices(); i++) {
        const Vertex& ve = P.vertex(i);
        // if (v == 761) {
        //     logger().debug("seed {} has polygon with vertex: ({},{},{})",
        //         ve.point()[0],ve.point()[1], ve.point()[2]
        //     );
        // }
        builder_.add_vertex_to_facet(ve.point(), ve.sym());
      }
      builder_.end_facet();

      GEO::Process::release_spinlock(global_lock_);

      // std::cout << "builder finish to build facet" << std::endl;
    }

   private:
    Process::spinlock global_lock_;
    const GenRestrictedPowerDiagram& RPD;
    BUILDER& builder_;
    signed_index_t current_facet_;
  };

  // Empty class doing nothing but print
  // this is for dubug only
  class DebugBuildRPD {
   public:
    /**
     * \brief Constructs a new BuildRPD.
     * \param[in] RPD_in the restricted Voronoi diagram
     * \param[in] builder the lesh builder
     */
    DebugBuildRPD() {
      std::cout << "initializing DebugBuildRPD ..." << std::endl;
    }

    /**
     * \brief The destructor
     * \details Terminates the current facet
     *    and the current surface.
     */
    ~DebugBuildRPD() {
      std::cout << "calling ~DebugBuildRPD() ..." << std::endl;
      std::cout << "finish ~DebugBuildRPD()" << std::endl;
    }

    /**
     * \brief The callback called for each restricted Voronoi cell.
     * \param[in] v index of current center vertex
     * \param[in] P current restricted Voronoi cell
     */
    void operator()(const GEO::index_t v,
                    const typename GenRestrictedPowerDiagram::Polygon& P) {
      std::cout << "In DebugBuildRPD operator ..." << std::endl;
      std::cout << "builder finish to build facet" << std::endl;
    }
  };

  // NOT USE!!!
  /**
   * \brief Implementation class for explicitly constructing
   *    a volumetric mesh that corresponds to the volumetric
   *    restricted Voronoi diagram.
   * \details To be used as a template argument
   *    to RPD::for_each_volumetric_integration_simplex().
   *    The current Voronoi cell is reported in tetrahedron region.
   * \note For the moment, vertices are duplicated (will be fixed
   *  in a future version).
   */
  class BuildVolumetricRPD {
   public:
    /**
     * Constructs a new BuildVolumetricRPD.
     * \param[in] RPD the volumetric restricted Voronoi diagram
     * \param[in] dim dimension of the points (can be smaller
     *  than actual dimension of the RPD).
     * \param[out] vertices coordinates of the generated vertices
     * \param[out] triangle_vertex_indices generated triangles, as
     *  vertex indices triplets
     * \param[out] tet_vertex_indices generated tetrahedra, as
     *  vertex indices 4-uples
     * \param[out] triangle_regions each generated triangle has
     *  a region index, that corresponds to the index of the
     *  Voronoi cell the triangle belongs to
     * \param[out] tet_regions each generated tetrahedron has
     *  a region index, that corresponds to the index of the
     *  Voronoi cell the tetrahedron belongs to
     * \param[in] cell_borders_only if true, only the surfacic
     *  borders of the volumetric cells are saved in the mesh, else
     *  volumetric cells are tetrahedralized.
     * \pre dim <= delaunay->dimension()
     */
    BuildVolumetricRPD(
        GenRestrictedPowerDiagram& RPD,
        std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_seed_adj,
        coord_index_t dim, vector<double>& vertices,
        vector<index_t>& triangle_vertex_indices,
        vector<index_t>& tet_vertex_indices,
        vector<signed_index_t>& triangle_regions,
        vector<signed_index_t>& triangle_adj_regions,
        vector<signed_index_t>& tet_regions, bool cell_borders_only)
        : rt_(RPD.delaunay()),
          mesh_(RPD.mesh()),
          rpd_point_adj_(rpd_seed_adj),
          dim_(dim),
          vertices_(vertices),
          triangle_vertex_indices_(triangle_vertex_indices),
          tet_vertex_indices_(tet_vertex_indices),
          triangle_regions_(triangle_regions),
          triangle_adj_regions_(triangle_adj_regions),
          tet_regions_(tet_regions),
          cell_borders_only_(cell_borders_only) {
      vertices_.clear();
      triangle_vertex_indices_.clear();
      tet_vertex_indices_.clear();
      triangle_regions_.clear();
      triangle_adj_regions_.clear();
      tet_regions_.clear();
      global_lock_ = GEOGRAM_SPINLOCK_INIT;

      // // The first vertices are copied from Delaunay,
      // // the other ones will be created during the traversal
      // nb_vertices_ = rt_->get_nb_vertices();
      // vertices_.resize(nb_vertices_ * dim);
      // for(GEO::index_t v = 0; v < nb_vertices_; ++v) {
      //     double* p = rt_->get_double_data(v);
      //     for(coord_index_t c = 0; c < dim; ++c) {
      //         vertices_[v * dim + c] = p[c];
      //     }
      // }

      // // NOTE: looping through rt_->get_nb_vertices() would cause error
      // //       when calling rt_->get_double_data(v)
      // //       but since all points are not surface points, we do not need
      // //       to save them
      nb_vertices_ = 0;
      vertex_map_.set_first_vertex_index(nb_vertices_);
    }

    /**
     * \brief The callback called for each integration simplex.
     * \param[in] v index of current center vertex
     * \param[in] v_adj (unused here) is the index of the Voronoi cell
     *  adjacent to t accros facet (\p v1, \p v2, \p v3) or
     *  -1 if it does not exists
     *  \param[in] t (unused here) is the index of the current
     *   tetrahedron
     *  \param[in] t_adj (unused here) is the index of the
     *   tetrahedron adjacent to t accros facet (\p v1, \p v2, \p v3)
     *   or -1 if it does not exists
     * \param[in] v1 first vertex of current integration simplex
     * \param[in] v2 second vertex of current integration simplex
     * \param[in] v3 third vertex of current integration simplex
     */
    void operator()(index_t v, signed_index_t v_adj, index_t t,
                    signed_index_t t_adj, const Vertex& v1, const Vertex& v2,
                    const Vertex& v3) {
      geo_argused(v_adj);
      geo_argused(t);

      // logger().debug("In integration simplex: seed {} has adjacent seed {}",
      // v, v_adj);

      index_t iv2 = find_or_create_vertex(v, v1);
      index_t iv3 = find_or_create_vertex(v, v2);
      index_t iv4 = find_or_create_vertex(v, v3);

      if (cell_borders_only_) {
        // if(signed_index_t(v) > v_adj) {
        triangle_vertex_indices_.push_back(iv4);
        triangle_vertex_indices_.push_back(iv3);
        triangle_vertex_indices_.push_back(iv2);
        triangle_regions_.push_back(signed_index_t(v));
        triangle_adj_regions_.push_back(v_adj);
        // }
      } else {
        // note:
        // tet contains iv1 = v, which means
        // we need to store v first in vertices_
        index_t iv1 = v;
        // index_t iv2 = find_or_create_vertex(v, v1);
        // index_t iv3 = find_or_create_vertex(v, v2);
        // index_t iv4 = find_or_create_vertex(v, v3);

        // Triangle v1,v2,v3 is on border if there is
        // no adjacent seed and no adjacent tet.
        if (v_adj == -1 && t_adj == -1) {
          triangle_vertex_indices_.push_back(iv4);
          triangle_vertex_indices_.push_back(iv3);
          triangle_vertex_indices_.push_back(iv2);
          triangle_regions_.push_back(signed_index_t(v));
        }

        tet_vertex_indices_.push_back(iv1);
        tet_vertex_indices_.push_back(iv2);
        tet_vertex_indices_.push_back(iv3);
        tet_vertex_indices_.push_back(iv4);
        tet_regions_.push_back(signed_index_t(v));
      }

      index_t f_idx = triangle_regions_.size() - 1;
      if (f_idx == 85612 || f_idx == 85611 || f_idx == 8260) {
        logger().debug("fidx {} has vertices ({},{},{})", f_idx, iv4, iv3, iv2);
      }

      // store adjacent info
      if (rpd_point_adj_ != nullptr && v_adj != -1 && v != v_adj) {
        const double* v1_p = v1.point();
        const double* v2_p = v2.point();
        const double* v3_p = v3.point();
        Vector3 v1_tmp(v1_p[0], v1_p[1], v1_p[2]);
        Vector3 v2_tmp(v2_p[0], v2_p[1], v2_p[2]);
        Vector3 v3_tmp(v3_p[0], v3_p[1], v3_p[2]);
        double area = get_triangle_area(v1_tmp, v2_tmp, v3_tmp);
        if (area < 0.01) {
          // // check if v_adj is stored, if so, delete
          // std::set<GEO::index_t>& adj_points = (*rpd_point_adj_)[v];
          // auto it = adj_points.find(v_adj);
          // if (it != adj_points.end()) {
          //     adj_points.erase(it);
          //     // logger().debug("tag {} deleted v_adj {}", v, v_adj);
          // }

          // logger().debug("face has area: {}, v {}, v_adj {}, t {}", area, v,
          // v_adj, t); return;
        }

        // Vector3 vec12 = v2_tmp - v1_tmp;
        // Vector3 vec13 = v3_tmp - v1_tmp;
        // double angle = angle_between_two_vectors_in_degrees(vec12, vec13);
        // if (angle < 3 || angle > 177) {
        //     // logger().debug("face has angle: {}, v {}, v_adj {}, t {}",
        //     angle, v, v_adj, t); return;
        // }

        // std::set<int> adj_tag_debug = {669, 1208, 1464, 1477};
        // if (v == 1456 && adj_tag_debug.find(v_adj) != adj_tag_debug.end()) {
        //     logger().debug("tag {} v_adj {} has face {} with area {}, t {}",
        //         v, v_adj, triangle_regions_.size()-1,
        //         area, t
        //     );
        //     logger().debug("face {} has vertices ({},{},{}), ({},{},{}),
        //     ({},{},{})",
        //         triangle_regions_.size()-1,
        //         v1_tmp[0], v1_tmp[1], v1_tmp[2],
        //         v2_tmp[0], v2_tmp[1], v2_tmp[2],
        //         v3_tmp[0], v3_tmp[1], v3_tmp[2]
        //     );
        // }

        GEO::Process::acquire_spinlock(global_lock_);
        (*rpd_point_adj_)[v].insert(v_adj);
        GEO::Process::release_spinlock(global_lock_);
      }
    }

    /**
     * \brief The callback called for each tetrahedron
     * \param[in] v index of current center vertex
     * \param[in] v_adj is the index of the Voronoi cell
     *  adjacent to t accros facet (\p v2, \p v3, \p v4) or
     *  -1 if it does not exists
     *  \param[in] t (unused here) is the index of the current
     *   tetrahedron
     *  \param[in] t_adj (unused here) is the index of the
     *   tetrahedron adjacent to t accros facet (\p v2, \p v3, \p v4)
     *   or -1 if it does not exists
     * \param[in] v1 first vertex of current tetrahedron
     * \param[in] v2 second vertex of current tetrahedron
     * \param[in] v3 third vertex of current tetrahedron
     * \param[in] v4 fourth vertex of current tetrahedron
     */
    void operator()(index_t v, signed_index_t v_adj, index_t t,
                    signed_index_t t_adj, const Vertex& v1, const Vertex& v2,
                    const Vertex& v3, const Vertex& v4) {
      geo_argused(v_adj);
      geo_argused(t);
      geo_argused(t_adj);

      // logger().debug("seed {} has adjacent seed {}", v, v_adj);

      index_t iv1 = vertices_.size() / dim_;
      for (index_t c = 0; c < dim_; ++c) {
        vertices_.push_back(v1.point()[c]);
      }
      index_t iv2 = vertices_.size() / dim_;
      for (index_t c = 0; c < dim_; ++c) {
        vertices_.push_back(v2.point()[c]);
      }
      index_t iv3 = vertices_.size() / dim_;
      for (index_t c = 0; c < dim_; ++c) {
        vertices_.push_back(v3.point()[c]);
      }
      index_t iv4 = vertices_.size() / dim_;
      for (index_t c = 0; c < dim_; ++c) {
        vertices_.push_back(v4.point()[c]);
      }
      tet_vertex_indices_.push_back(iv1);
      tet_vertex_indices_.push_back(iv2);
      tet_vertex_indices_.push_back(iv3);
      tet_vertex_indices_.push_back(iv4);
      tet_regions_.push_back(signed_index_t(v));

      // NOTE: somehow this cannot return us the proper
      //       neighboring info
      if (rpd_point_adj_ != nullptr && v_adj != -1 && v != v_adj) {
        GEO::Process::acquire_spinlock(global_lock_);
        (*rpd_point_adj_)[v].insert(v_adj);
        GEO::Process::release_spinlock(global_lock_);
      }
    }

   protected:
    /**
     * \brief Retrieves the index of a vertex given its symbolic
     * representation.
     * \param[in] center_vertex_id index of current Voronoi seed
     * \param[in] v symbolic and geometric representation of the vertex
     * \return the index of the vertex
     */
    index_t find_or_create_vertex(index_t center_vertex_id, const Vertex& v) {
      index_t result =
          vertex_map_.find_or_create_vertex(center_vertex_id, v.sym());
      if (result == 231992 || result == 231989 || result == 40121) {
        logger().debug("new vertex {} is created, tag {}", result,
                       center_vertex_id);
      }

      if (result >= nb_vertices_) {
        geo_assert(result == nb_vertices_);
        nb_vertices_ = result + 1;
        for (coord_index_t c = 0; c < dim_; ++c) {
          vertices_.push_back(v.point()[c]);
        }
      }
      return result;
    }

   private:
    // global lock for updating rpd_point_adj_
    Process::spinlock global_lock_;
    std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_point_adj_;

    RegularTriangulationNN* rt_;
    const GEO::Mesh* mesh_;
    coord_index_t dim_;
    vector<double>& vertices_;
    vector<index_t>& triangle_vertex_indices_;
    vector<index_t>& tet_vertex_indices_;
    vector<signed_index_t>& triangle_regions_;
    vector<signed_index_t>& triangle_adj_regions_;
    vector<signed_index_t>& tet_regions_;
    matfp::RPDVertexMap vertex_map_;
    index_t nb_vertices_;
    bool cell_borders_only_;
  };  // end BuildVolumetricRPD

  class BuildRPDVolumetricMesh {
   private:
    // facet degenerated type
    enum DEG_TYPE {
      NONE = 0, /*not degenerate*/
      POINT = 1,
      LINE_2_1 = 2,
      LINE_1_1_1 = 3
    };

   public:
    /**
     * Constructs a new BuildVolumetricRPD.
     * \param[in] RPD the volumetric restricted Voronoi diagram
     * \param[in] dim dimension of the points (can be smaller
     *  than actual dimension of the RPD).
     * \param[out] vertices coordinates of the generated vertices
     * \param[out] triangle_vertex_indices generated triangles, as
     *  vertex indices triplets
     * \param[out] tet_vertex_indices generated tetrahedra, as
     *  vertex indices 4-uples
     * \param[out] triangle_regions each generated triangle has
     *  a region index, that corresponds to the index of the
     *  Voronoi cell the triangle belongs to
     * \param[out] tet_regions each generated tetrahedron has
     *  a region index, that corresponds to the index of the
     *  Voronoi cell the tetrahedron belongs to
     * \param[in] cell_borders_only if true, only the surfacic
     *  borders of the volumetric cells are saved in the mesh, else
     *  volumetric cells are tetrahedralized.
     * \pre dim <= delaunay->dimension()
     */
    BuildRPDVolumetricMesh(
        GenRestrictedPowerDiagram& RPD,
        std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_seed_adj,
        coord_index_t dim, vector<double>& vertices,
        vector<index_t>& triangle_vertex_indices,
        vector<index_t>& tet_vertex_indices,
        vector<signed_index_t>& triangle_regions,
        vector<signed_index_t>& triangle_adj_regions,
        vector<signed_index_t>& tet_regions, bool cell_borders_only)
        : rt_(RPD.delaunay()),
          mesh_(RPD.mesh()),
          rpd_point_adj_(rpd_seed_adj),
          dim_(dim),
          vertices_(vertices),
          triangle_vertex_indices_(triangle_vertex_indices),
          tet_vertex_indices_(tet_vertex_indices),
          triangle_regions_(triangle_regions),
          triangle_adj_regions_(triangle_adj_regions),
          tet_regions_(tet_regions),
          cell_borders_only_(cell_borders_only) {
      vertices_.clear();
      triangle_vertex_indices_.clear();
      tet_vertex_indices_.clear();
      triangle_regions_.clear();
      triangle_adj_regions_.clear();
      tet_regions_.clear();
      global_lock_ = GEOGRAM_SPINLOCK_INIT;

      rpd_point_adj_->clear();

      // Note: since all mats are not surface points, we do not need
      //       to save them
      nb_vertices_ = 0;
      vertex_map_.set_first_vertex_index(nb_vertices_);
    }

    /**
     * \brief The callback called for each integration simplex.
     * \param[in] v index of current center vertex
     * \param[in] v_adj (unused here) is the index of the Voronoi cell
     *  adjacent to t accros facet (\p v1, \p v2, \p v3) or
     *  -1 if it does not exists
     *  \param[in] t (unused here) is the index of the current
     *   tetrahedron
     *  \param[in] t_adj (unused here) is the index of the
     *   tetrahedron adjacent to t accros facet (\p v1, \p v2, \p v3)
     *   or -1 if it does not exists
     * \param[in] v1 first vertex of current integration simplex
     * \param[in] v2 second vertex of current integration simplex
     * \param[in] v3 third vertex of current integration simplex
     */
    void operator()(index_t v, signed_index_t v_adj, index_t t,
                    signed_index_t t_adj, const Vertex& v1, const Vertex& v2,
                    const Vertex& v3) {
      geo_argused(v_adj);
      geo_argused(t);

      // if(signed_index_t(v) <= v_adj) {
      //     return;
      // }

      // logger().debug("In integration simplex: seed {} has adjacent seed {}",
      // v, v_adj); logger().debug("processing seed {} and adj seed {}", v,
      // v_adj);

      index_t iv2 = find_or_create_vertex(v, v1);
      index_t iv3 = find_or_create_vertex(v, v2);
      index_t iv4 = find_or_create_vertex(v, v3);

      std::array<index_t, 3> vs_indices = {{iv2, iv3, iv4}};
      std::sort(vs_indices.begin(), vs_indices.end());

      // global lock
      GEO::Process::acquire_spinlock(global_lock_);

      // only store facet with type DEG_TYPE::NONE
      // triangle_vertex_indices_.push_back(iv4);
      // triangle_vertex_indices_.push_back(iv3);
      // triangle_vertex_indices_.push_back(iv2);
      triangle_vertex_indices_.push_back(iv2);
      triangle_vertex_indices_.push_back(iv3);
      triangle_vertex_indices_.push_back(iv4);
      triangle_regions_.push_back(signed_index_t(v));
      triangle_adj_regions_.push_back(signed_index_t(v_adj));

      // index_t f_idx = triangle_regions_.size()-1;
      // if (f_idx == 85612 || f_idx == 85611 || f_idx == 8260) {
      //     logger().debug("fidx {} has vertices ({},{},{})", f_idx, iv4, iv3,
      //     iv2);
      // }

      // store adjacent info
      if (rpd_point_adj_ != nullptr && v_adj != -1 && v != v_adj) {
        (*rpd_point_adj_)[v].insert(v_adj);
      }
      GEO::Process::release_spinlock(global_lock_);
    }

   protected:
    /**
     * \brief Retrieves the index of a vertex given its symbolic
     * representation.
     * \param[in] center_vertex_id index of current Voronoi seed
     * \param[in] v symbolic and geometric representation of the vertex
     * \return the index of the vertex
     */
    index_t find_or_create_vertex(index_t center_vertex_id, const Vertex& v) {
      index_t result =
          vertex_map_.find_or_create_vertex(center_vertex_id, v.sym());

      if (result >= nb_vertices_) {
        geo_assert(result == nb_vertices_);
        nb_vertices_ = result + 1;
        for (coord_index_t c = 0; c < dim_; ++c) {
          vertices_.push_back(v.point()[c]);
        }
      }
      return result;
    }

   private:
    RegularTriangulationNN* rt_;
    const GEO::Mesh* mesh_;
    coord_index_t dim_;
    vector<double>& vertices_;
    vector<index_t>& triangle_vertex_indices_;
    vector<index_t>& tet_vertex_indices_;
    vector<signed_index_t>& triangle_regions_;
    vector<signed_index_t>& triangle_adj_regions_;
    vector<signed_index_t>& tet_regions_;
    matfp::RPDVertexMap vertex_map_;
    index_t nb_vertices_;
    bool cell_borders_only_;

    // global lock for updating rpd_point_adj_
    GEO::Process::spinlock global_lock_;
    std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_point_adj_;
  };

  /********************************************************************/

  /**
   * Helper function for building RPD surfacic mesh in parallel
   **/
  virtual void build_rpd_mesh_surfacic(BuildRPD<RPDMeshBuilder>& build_rpd) {
    // logger().debug("calling build_rpd_mesh_surfacic ...");
    create_threads();
    if (nb_parts() == 0) {
      // logger().debug("in build_rpd_mesh_surfacic nb_parts() == 0");
      // PolygonCallbackAction action(RPD_,polygon_callback);
      RPD_.for_each_polygon(build_rpd);
    } else {
      // logger().debug("in build_rpd_mesh_surfacic nb_parts() != 0");

      for (index_t t = 0; t < nb_parts(); t++) {
        part(t).RPD_.set_symbolic(RPD_.symbolic());
        part(t).RPD_.set_connected_components_priority(
            RPD_.connected_components_priority());
      }
      // Note: we need global lock here, managed by BuildRPD
      thread_mode_ = MT_RPD_S_MESH;
      build_rpd_ = &build_rpd;
      // Note: callback begin()/end() is called in for_each_polygon()
      parallel_for(0, nb_parts(), [this](index_t i) { run_thread(i); });
    }
  }

  void compute_RPD(
      GEO::Mesh& M,
      std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_seed_adj,
      std::map<GEO::index_t, std::set<GEO::index_t>>* rpd_vs_bisectors,
      GEO::coord_index_t dim, bool cell_borders_only,
      bool integration_simplices, bool is_parallel) override {
    bool sym = RPD_.symbolic();
    RPD_.set_symbolic(true);

    if (volumetric_) {
      logger().debug("ERROR: cannot compute volumetric RPD.");
      log_and_throw("ERROR");
    } else {  // surfacic

      matfp::RPDMeshBuilder builder(&M, mesh_, rpd_seed_adj, rpd_vs_bisectors);
      if (dim != 0) {
        builder.set_dimension(dim);
      }
      BuildRPD<RPDMeshBuilder> build_rpd_action(RPD_, builder);

      if (is_parallel) {
        logger().debug("computing RPD surfacic in parallel ...");
        build_rpd_mesh_surfacic(build_rpd_action);
      } else {
        logger().debug("computing RPD surfacic in sequence ...");
        RPD_.for_each_polygon(build_rpd_action);
      }

      // RPD_.for_each_polygon(
      //     DebugBuildRPD()
      // );

      // logger().debug("after for_each_polygon");
      // std::cout << "in compute_RPD, finish for_each_polygon" << std::endl;
    }

    RPD_.set_symbolic(sym);
    M.show_stats("RPD");
  }

  /********************************************************************/
  /**
   * \brief Place holder, "no locking" policy.
   * \details NoLocks is used by algorithms templated
   *  by locking policy, for the single-threaded instances
   *  that do not need synchronization. The multi-threaded
   *  instances are parameterized by SpinLockArray.
   */
  class NoLocks {
   public:
    /**
     * \brief Acquires a spinlock.
     * \details Does nothing in this version
     * \param[in] i index of the spinlock to acquire
     */
    void acquire_spinlock(index_t i) { GEO::geo_argused(i); }

    /**
     * \brief Releases a spinlock.
     * \details Does nothing in this version
     * \param[in] i index of the spinlock to release
     */
    void release_spinlock(index_t i) { GEO::geo_argused(i); }
  };

  /********************************************************************/

  void for_each_polygon(matfp::RPDPolygonCGALCallback& callback, bool symbolic,
                        bool connected_comp_priority, bool parallel) override {
    bool sym_backup = RPD_.symbolic();
    RPD_.set_symbolic(symbolic);
    RPD_.set_connected_components_priority(connected_comp_priority);
    callback.begin();
    if (parallel) {
      compute_with_polygon_callback(callback);
    } else {
      PolygonCallbackAction action(RPD_, callback);
      RPD_.for_each_polygon(action);
    }
    callback.end();
    RPD_.set_symbolic(sym_backup);
    RPD_.set_connected_components_priority(false);
  }

  void for_each_polyhedron(matfp::RPDPolyhedronCGALCallback& callback,
                           bool symbolic, bool connected_comp_priority,
                           bool parallel) override {
    bool sym_backup = RPD_.symbolic();
    RPD_.set_symbolic(symbolic);
    RPD_.set_connected_components_priority(connected_comp_priority);
    callback.set_dimension(RPD_.mesh()->vertices.dimension());
    callback.begin();
    if (parallel) {
      compute_with_polyhedron_callback(callback);
    } else {
      RPD_.for_each_polyhedron(callback);
    }
    callback.end();
    RPD_.set_symbolic(sym_backup);
    RPD_.set_connected_components_priority(false);
  }

  /********************************************************************/
  /**
   * \brief Does the actual computation for a specific part
   *    in multithread mode.
   * \param[in] t the index of the part.
   * \pre \p t < nb_parts()
   */
  void run_thread(index_t t) {
    geo_assert(t < nb_parts());
    thisclass& T = part(t);
    switch (thread_mode_) {
      // case MT_LLOYD:
      // {
      //     T.compute_centroids(arg_vectors_, arg_scalars_);
      // } break;
      // case MT_NEWTON:
      // {
      //     T.compute_CVT_func_grad(T.funcval_, arg_vectors_);
      // } break;
      // case MT_INT_SMPLX:
      // {
      //     T.compute_integration_simplex_func_grad(
      //         T.funcval_, arg_vectors_, simplex_func_
      //     );
      // } break;
      case MT_POLYG: {
        T.compute_with_polygon_callback(*polygon_callback_);
      } break;
      case MT_POLYH: {
        T.compute_with_polyhedron_callback(*polyhedron_callback_);
      } break;
      case MT_RPD_S_MESH: {
        T.build_rpd_mesh_surfacic(*build_rpd_);
      } break;
      case MT_NONE:
        geo_assert_not_reached;
    }
  }

  /********************************************************************/

  void set_delaunay(RegularTriangulationNN* rt) override {
    baseclass::set_delaunay(rt);
    RPD_.set_delaunay(rt);
    for (index_t p = 0; p < nb_parts_; ++p) {
      parts_[p].set_delaunay(rt);
    }
  }

  void set_check_SR(bool x) override {
    RPD_.set_check_SR(x);
    for (index_t p = 0; p < nb_parts_; ++p) {
      parts_[p].set_check_SR(x);
    }
  }

  void set_exact_predicates(bool x) override {
    RPD_.set_exact_predicates(x);
    for (index_t p = 0; p < nb_parts_; ++p) {
      parts_[p].set_exact_predicates(x);
    }
  }

  void create_threads() override {
    // TODO: check if number of facets is not smaller than
    // number of threads
    // TODO: create parts even if facets range is specified
    // (and subdivide facets range)
    if (is_slave_ || facets_begin_ != -1 || facets_end_ != -1) {
      return;
    }
    index_t nb_parts_in = Process::maximum_concurrent_threads();
    if (nb_parts() != nb_parts_in) {
      if (nb_parts_in == 1) {
        delete_threads();
      } else {
        vector<index_t> facet_ptr;
        vector<index_t> tet_ptr;
        mesh_partition(*mesh_, MESH_PARTITION_HILBERT, facet_ptr, tet_ptr,
                       nb_parts_in);
        delete_threads();
        parts_ = new thisclass[nb_parts_in];
        nb_parts_ = nb_parts_in;
        for (index_t i = 0; i < nb_parts(); ++i) {
          part(i).mesh_ = mesh_;
          part(i).set_delaunay(rt_);
          part(i).R3_embedding_base_ = R3_embedding_base_;
          part(i).R3_embedding_stride_ = R3_embedding_stride_;
          part(i).master_ = this;
          part(i).RPD_.set_mesh(mesh_);
          part(i).set_facets_range(facet_ptr[i], facet_ptr[i + 1]);
          part(i).set_exact_predicates(RPD_.exact_predicates());
          // part(i).set_volumetric(volumetric());
          part(i).set_check_SR(RPD_.check_SR());
        }
        // if(mesh_->cells.nb() != 0) {
        //     for(index_t i = 0; i < nb_parts(); ++i) {
        //         part(i).set_tetrahedra_range(
        //             tet_ptr[i], tet_ptr[i + 1]
        //         );
        //     }
        // }
        geo_assert(!Process::is_running_threads());
      }
    }
  }

  void set_volumetric(bool x) override {
    volumetric_ = x;
    for (index_t i = 0; i < nb_parts(); ++i) {
      part(i).set_volumetric(x);
    }
  }

  void set_facets_range(index_t facets_begin, index_t facets_end) override {
    RPD_.set_facets_range(facets_begin, facets_end);
    facets_begin_ = signed_index_t(facets_begin);
    facets_end_ = signed_index_t(facets_end);
  }

  void delete_threads() override {
    delete[] parts_;
    parts_ = nullptr;
    nb_parts_ = 0;
  }

  /**
   * \brief Gets the number of parts (or number of threads).
   */
  index_t nb_parts() const { return nb_parts_; }

  /**
   * \brief Gets a given part from its index.
   * \param[in] i index of the part
   * \pre \p i < nb_parts()
   */
  thisclass& part(index_t i) {
    geo_debug_assert(i < nb_parts());
    return parts_[i];
  }

  /**
   * \copydoc RestrictedVoronoiDiagram::point_allocator()
   */
  GEOGen::PointAllocator* point_allocator() override {
    return RPD_.point_allocator();
  }

 protected:
  GenRestrictedPowerDiagram RPD_;

  // For projection
  bool use_exact_projection_;
  index_t nb_triangles_;
  vector<index_t> triangles_;
  vector<vector<index_t>> stars_;
  Delaunay_var mesh_vertices_;

  // One of MT_NONE, MT_LLOYD, MT_NEWTON
  ThreadMode thread_mode_;

  bool is_slave_;

  // Variables for 'master' in multithreading mode
  thisclass* parts_;
  index_t nb_parts_;
  Process::SpinLockArray spinlocks_;

  // Newton mode with int. simplex
  IntegrationSimplex* simplex_func_;

  // PolygonCallback mode.
  RPDPolygonCGALCallback* polygon_callback_;

  // PolyhedronCallback mode.
  // NO implementation, use RVD
  RPDPolyhedronCGALCallback* polyhedron_callback_;

  BuildRPD<RPDMeshBuilder>* build_rpd_;

  // master stores argument for compute_centroids() and
  // compute_CVT_func_grad() to pass it to the parts.
  double* arg_vectors_;
  double* arg_scalars_;

  // Variables for 'slaves' in multithreading mode
  thisclass* master_;
  double funcval_;  // Newton mode: function value

 protected:
  /**
   * \brief Destructor
   */
  ~RPD_3d_Impl() override { delete_threads(); }

 private:
  /** \brief Forbids construction by copy. */
  RPD_3d_Impl(const thisclass&);

  /** \brief Forbids assignment. */
  thisclass& operator=(const thisclass&);
};

}  // namespace

namespace matfp {

RestrictedPowerDiagram* RestrictedPowerDiagram::create(
    RegularTriangulationNN* rt, GEO::Mesh* mesh, const double* R3_embedding,
    index_t R3_embedding_stride) {
  geo_cite("DBLP:journals/tog/EdelsbrunnerM90");
  geo_cite("DBLP:conf/compgeom/Shewchuk96");
  geo_cite("meyer:inria-00344297");
  geo_cite("DBLP:conf/gmp/YanWLL10");
  geo_cite("DBLP:journals/cad/YanWLL13");
  geo_cite("DBLP:journals/cad/Levy16");

  // delaunay->set_stores_neighbors(true);
  RestrictedPowerDiagram* result = nullptr;
  geo_assert(rt != nullptr);

  result = new RPD_3d_Impl(rt, mesh, R3_embedding, R3_embedding_stride);

  // if(CmdLine::get_arg("algo:predicates") == "exact") {
  result->set_exact_predicates(true);
  // }
  return result;
}

RestrictedPowerDiagram::~RestrictedPowerDiagram() {}

RestrictedPowerDiagram::RestrictedPowerDiagram(RegularTriangulationNN* rt,
                                               Mesh* mesh,
                                               const double* R3_embedding,
                                               index_t R3_embedding_stride)
    : dimension_(0),
      mesh_(mesh),
      R3_embedding_base_(R3_embedding),
      R3_embedding_stride_(R3_embedding_stride) {
  set_delaunay(rt);
  // has_weights_ = false;
  facets_begin_ = -1;
  facets_end_ = -1;
  tets_begin_ = -1;
  tets_end_ = -1;
  // volumetric_ = false;
}

void RestrictedPowerDiagram::set_delaunay(RegularTriangulationNN* rt) {
  rt_ = rt;
  if (rt_ != nullptr) {
    dimension_ = 3;
  } else {
    dimension_ = 0;
  }
}

}  // namespace matfp