// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <matfp/Logger.h>
#include <matfp/Types/CommonTypes.h>

namespace matfp {

class Cone {
 public:
  enum Type {
    INVALID = 1,   // -- height = 0.0
    LINE = 2,      // -- line (equal radii & radius = 0.0)
    CYLINDER = 3,  // -- cylinder (equal radii)
    CONE = 4       // -- general regular cone (vary radii)
  };

 public:
  //////////
  // We define cone as:
  // 	base radius <= top radius

  // parameters used by Gmsh
  Vector3 bcenter = Vector3(0, 0, 1);     // center of base cirvular face
  Vector3 tcenter = Vector3(0, 0, 1);     // center of top cirvular face
  Vector3 diffcenter = Vector3(0, 0, 0);  // (dx, dy, dz): tcenter - bcenter
  double base;                            // radius of base circular face
  double top;                             // radius of top circular face

  // 1 -- height = 0.0
  // 2 -- line (equal radii & radius = 0.0)
  // 3 -- cylinder (equal radii)
  // 4 -- general regular cone (vary radii)
  Type type;

  // parameters used by gluCylinder()
  // we also need rotation axis&angle to mesh a cone by hand
  // check Envelope::generate_one_cone()
  Vector3 apex;  // todo: not using
  Vector3 axis;  // unit vector from base to top
  double height;
  // the rotation axis to rotate the axis to z-axis
  Vector3 rot_axis;
  double rot_angle;          // angle in degrees
  double rot_angle_radians;  // angle in radians

 public:
  Cone() {}
  Cone(Vector3 c0, double r0, Vector3 c1, double r1);

 public:
  inline void print_info() {
    logger().info("------ Cone Info ------");
    logger().info("type: {}", static_cast<int>(type));
    logger().info("bcenter: ({},{},{})", bcenter[0], bcenter[1], bcenter[2]);
    logger().info("tcenter: ({},{},{})", tcenter[0], tcenter[1], tcenter[2]);
    logger().info("diffcenter: ({},{},{})", diffcenter[0], diffcenter[1],
                  diffcenter[2]);
    logger().info("base: {}", base);
    logger().info("top: {}", top);
    logger().info("height: {}", height);
    logger().info("apex: ({},{},{})", apex[0], apex[1], apex[2]);
    logger().info("axis: ({},{},{})", axis[0], axis[1], axis[2]);
    logger().info("rot_angle: {}", rot_angle);
    logger().info("rot_axis: ({},{},{})", rot_axis[0], rot_axis[1],
                  rot_axis[2]);
  }
};

class SimpleTriangle {
 public:
  SimpleTriangle() : normal(Vector3::Zero()) {};

 public:
  inline void print_st_info() {
    logger().debug("-------- Simple Triangle Info -------");
    for (int i = 0; i < 3; i++)
      logger().debug("v[{}]: ({},{},{})", i, v[i][0], v[i][1], v[i][2]);
    logger().debug("normal: ({},{},{})", normal[0], normal[1], normal[2]);
  }

 public:
  ///////////
  // the orientaton of verices is counter-clokwise
  Vector3 v[3];
  Vector3 normal;

  // bool project_point_onto_triangle(
  //     const Vector3& p,
  //     Vector3& p_proj,
  //     double& dist
  // );
  void update_normal(bool is_reverse = false);
  //////////
  // This update would make sure all vertices are counter-clockwise
  // and normal must be pointing outside
  void update_normal(Vector3& _normal);
};

// return types:
// 1 -- INVALID
// 2 -- TETRAHEDRON
// 3 -- PRISM
// here matching NonManifoldMesh_Face::Type
int get_triangles_from_three_spheres(const Vector3& c0, const double& r0,
                                     const Vector3& c1, const double& r1,
                                     const Vector3& c2, const double& r2,
                                     SimpleTriangle& st0, SimpleTriangle& st1);

}  // namespace matfp
