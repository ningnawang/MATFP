// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include "Primitives.h"

#include <matfp/Common.h>

namespace matfp {

Cone::Cone(Vector3 c0, double r0, Vector3 c1, double r1) {
  Vector3 c0c1 = c1 - c0;
  // one sphere is included in another sphere
  if (c0c1.norm() - std::abs(r1 - r0) < SCALAR_ZERO) {
    apex = r1 > r0 ? c1 : c0;
    axis = Vector3(0, 0, 1);
    base = r1 > r0 ? r1 : r0;
    top = r1 > r0 ? r1 : r0;
    height = 0.;
    type = Cone::Type::INVALID;
    rot_axis = Vector3(0, 0, 1);
    rot_angle = 0.;
    rot_angle_radians = 0.;
    return;
  }
  if (c0c1.norm() < SCALAR_ZERO) {
    apex = c0;
    axis = Vector3(0, 0, 1);
    base = r0;
    top = r0;
    height = 0.;
    type = Cone::Type::INVALID;
    rot_axis = Vector3(0, 0, 1);
    rot_angle = 0.;
    rot_angle_radians = 0.;

  } else {
    double dr0r1 = std::fabs(r0 - r1);
    // cylinder
    if (dr0r1 < SCALAR_ZERO) {
      bcenter = c0;
      tcenter = c1;
      diffcenter = tcenter - bcenter;
      base = r0;
      top = r0;

      apex = c0;
      axis = c1 - c0;
      axis.normalize();
      height = (c1 - c0).norm();

      if (base < SCALAR_ZERO) {
        type = Cone::Type::LINE;
      } else {
        type = Cone::Type::CYLINDER;
      }
      // normal cones (vary radii)
    } else {
      apex = (r1 * c0 - r0 * c1) / (r1 - r0);
      axis = (r0 < r1) ? (c1 - c0) : (c0 - c1);
      axis.normalize();

      // logger().debug("-----------------------");
      // logger().debug("apex: ({},{},{})", apex[0], apex[1], apex[2]);
      // logger().debug("c0: ({},{},{})", c0[0], c0[1], c0[2]);
      // logger().debug("c1: ({},{},{})", c1[0], c1[1], c1[2]);
      // logger().debug("r0: {}", r0);
      // logger().debug("r1: {}", r1);

      double cangle;
      Vector3 apexc0 = c0 - apex;
      double vc0len = apexc0.norm();
      Vector3 apexc1 = c1 - apex;
      double vc1len = apexc1.norm();
      // cangle = std::sqrt(1.-r0*r0/vc0len/vc0len);
      // NOTE: r0 or r1 might be 0.0, but not both (Cone::Type::LINE if both)
      // 		 use the std::max(r0,r1) to compute cangle
      if (r0 < r1) {
        cangle = std::sqrt(1. - r1 * r1 / vc1len / vc1len);
        bcenter = apex + apexc0 * cangle * cangle;
        tcenter = apex + apexc1 * cangle * cangle;
        diffcenter = tcenter - bcenter;
        base = r0 * cangle;
        top = r1 * cangle;
        // logger().debug("cangle: {}", cangle);
        // logger().debug("bcenter: ({},{},{})", bcenter[0], bcenter[1], c1[2]);

        apex = bcenter;
        height = (vc1len - vc0len) * cangle * cangle;
      } else {
        cangle = std::sqrt(1. - r0 * r0 / vc0len / vc0len);
        bcenter = apex + apexc1 * cangle * cangle;
        tcenter = apex + apexc0 * cangle * cangle;
        diffcenter = tcenter - bcenter;
        base = r1 * cangle;
        top = r0 * cangle;

        apex = bcenter;
        height = (vc0len - vc1len) * cangle * cangle;
      }
      type = Cone::Type::CONE;
    }

    Vector3 za(0, 0, 1);
    rot_angle_radians = std::acos(axis.dot(za));
    if ((std::fabs(rot_angle_radians) < SCALAR_ZERO) ||
        (std::fabs(rot_angle_radians - PI) < SCALAR_ZERO))
      rot_axis = Vector3(1, 0, 0);
    else
      rot_axis = axis.cross(za);
    rot_axis.normalize();
    rot_angle = rot_angle_radians * (180. / PI);

    // // sanity check
    // Vector3 diffcenter_norm = diffcenter.normalized();
    // if (!is_vector_same(diffcenter_norm, axis)) {
    // 	print_info();
    // 	logger().debug("c0: ({},{},{}), r0: {}, c1: ({},{},{}), r1: {}",
    // 		c0[0], c0[0], c0[0], r0,
    // 		c1[0], c1[0], c1[0], r1
    // 	);
    // 	logger().debug("diffcenter: ({},{},{}), diffcenter_norm: ({},{},{}),
    // axis: ({},{},{})", 		diffcenter[0], diffcenter[1],
    // diffcenter[2], 		diffcenter_norm[0], diffcenter_norm[1],
    // diffcenter_norm[2], 		axis[0], axis[1], axis[2]
    // 	);
    // 	log_and_throw("ERROR: diffcenter_norm and axis should be the same!!!!");
    // }
  }
}

void SimpleTriangle::update_normal(bool is_reverse) {
  normal = get_normal(v[0], v[1], v[2]);
  if (is_reverse) normal *= -1.;

  // Vector3 v01, v02;
  // v01 = v[1] - v[0];
  // v02 = v[2] - v[0];
  // normal = v01.cross(v02);
  // normal.normalize();
  // if (is_reverse)
  //     normal *= -1.;
}

// vs must be counter-clockwise, normal pointing outside
void SimpleTriangle::update_normal(Vector3& normal_) {
  // normal = normal_;
  Vector3 normal__ = get_normal(v[0], v[1], v[2]);
  normal = normal__;
  // there might be a little bit of distortion
  // of normal_ and normal__, but not that big.
  // the update here is mainly for checking orientations
  if (normal_.dot(normal__) >= 0) return;
  // normal__ and normal_ are opposite,
  // then swap vertices
  // logger().error("------------- SimpleTriangle: swapping vertices");
  Vector3 tmp = v[1];
  v[1] = v[2];
  v[2] = tmp;
  normal *= -1.;
  return;

  // if ((normal__ -normal_).norm() <= SCALAR_ZERO)
  // 	return;
  // // normal__ and normal_ are opposite
  // // swap vertices
  // if ((normal__ + normal_).norm() <= SCALAR_ZERO) {
  // 	// logger().error("------------- SimpleTriangle: swapping vertices");
  // 	Vector3 tmp = v[1];
  // 	v[1] = v[2];
  // 	v[2] = tmp;
  // 	return;
  // } else {
  // 	logger().error("ERROR: st normal ({},{},{}) not matching given normal
  // ({},{},{})", 		normal__[0], normal__[1], normal__[2],
  // normal_[0], normal_[1], normal_[2]
  // 	);
  // 	log_and_throw("ERROR: SimpleTriangle update normal failed");
  // }
}

int get_triangles_from_three_spheres(const Vector3& c0, const double& r0,
                                     const Vector3& c1, const double& r1,
                                     const Vector3& c2, const double& r2,
                                     SimpleTriangle& st0, SimpleTriangle& st1) {
  int type = 3;  // default as prism
  Vector3 c0c1(c1 - c0), c0c2(c2 - c0), c1c2(c2 - c1);
  double c0c1len(c0c1.norm()), c0c2len(c0c2.norm()), c1c2len(c1c2.norm());
  double dr0r1(std::fabs(r0 - r1)), dr0r2(std::fabs(r0 - r2)),
      dr1r2(std::fabs(r1 - r2));

  // some spheres are concentric and there are no triangles.
  if ((c0c1len < SCALAR_ZERO) || (c0c2len < SCALAR_ZERO) ||
      (c1c2len < SCALAR_ZERO)) {
    type = 1;
    return type;
  }

  //// some spheres are included in some other spheres
  // if ((c0c1len - abs(r0 - r1) < SCALAR_ZERO) || (c0c2len - abs(r0 - r2) <
  // SCALAR_ZERO) || (c1c2len - abs(r1 - r2) < SCALAR_ZERO)) 	return false;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ATTENTION: Orientation MATTERS!!
  // We obtain two simple triangles by moving c0,c1,c2 toward their normals.
  // Therefore, if we do NOT update orientatons, st0.v[i] and st1.v[i] must be
  // on the opposite side of a same sphere, i = 1,2,3 Since we update the order
  // of st0.v and st1.v by their normals to make sure orientation is alaways
  // counter-clockwise, the new matching pairs changed to: st0.v[0] - st1.v[0]
  // st0.v[1] - st1.v[2]
  // st0.v[2] - st1.v[1]
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  Vector3 norm;
  norm = c0c1.cross(c0c2);
  norm.normalize();

  // equal-radius spheres
  if ((dr0r1 < SCALAR_ZERO) && (dr0r2 < SCALAR_ZERO) && (dr1r2 < SCALAR_ZERO)) {
    st0.v[0] = c0 + norm * r0;
    st0.v[1] = c1 + norm * r1;
    st0.v[2] = c2 + norm * r2;
    st0.update_normal(norm);

    st1.v[0] = c0 - norm * r0;
    st1.v[1] = c1 - norm * r1;
    st1.v[2] = c2 - norm * r2;
    Vector3 norm1 = norm * -1.;
    st1.update_normal(norm1);
    return type;
  } else {
    // two points on the tangent plane
    Vector3 apex0, apex1;

    // two spheres are equal-radius
    if (dr0r1 < SCALAR_ZERO) {
      apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
      apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
      if (r0 < SCALAR_ZERO) type = 2;
    } else if (dr0r2 < SCALAR_ZERO) {
      apex0 = (r1 * c0 - r0 * c1) / (r1 - r0);
      apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
      if (r0 < SCALAR_ZERO) type = 2;
    } else if (dr1r2 < SCALAR_ZERO) {
      apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
      apex1 = (r0 * c1 - r1 * c0) / (r0 - r1);
      if (r1 < SCALAR_ZERO) type = 2;
    } else {
      apex0 = (r2 * c0 - r0 * c2) / (r2 - r0);
      apex1 = (r2 * c1 - r1 * c2) / (r2 - r1);
    }

    double distc0;
    Vector3 fp;
    double r_temp = r0;
    Vector3 c_temp;

    if ((apex0 - apex1).norm() <= SCALAR_ZERO) {
      if (dr0r1 < SCALAR_ZERO) {
        apex1 = c0c1 + apex0;
        project_point_onto_line(c1, apex0, apex1, fp, distc0);
        r_temp = r1;
        c_temp = c1;
      } else if (dr0r2 < SCALAR_ZERO) {
        apex1 = c0c2 + apex0;
        project_point_onto_line(c0, apex0, apex1, fp, distc0);
        r_temp = r0;
        c_temp = c0;
      } else {
        apex1 = c1c2 + apex0;
        project_point_onto_line(c2, apex0, apex1, fp, distc0);
        r_temp = r2;
        c_temp = c2;
      }
    } else if (r0 == 0) {
      if (r1 != 0) {
        project_point_onto_line(c1, apex0, apex1, fp, distc0);
        r_temp = r1;
        c_temp = c1;
      } else {
        project_point_onto_line(c2, apex0, apex1, fp, distc0);
        r_temp = r2;
        c_temp = c2;
      }
    } else {
      project_point_onto_line(c0, apex0, apex1, fp, distc0);
      r_temp = r0;
      c_temp = c0;
    }

    double sangle = r_temp / distc0;
    if (fabs(sangle) > 1.) {
      // std::cout<<"fabs(sangle): "<<fabs(sangle)<<std::endl;
      // std::cout<<"apex: "<<apex0<<" "<<apex1<<std::endl;
      type = 1;
      return type;
    }

    double cangle = sqrt(1. - r_temp * r_temp / distc0 / distc0);
    Vector3 norfpc0(c_temp - fp);
    norfpc0.normalize();

    Vector3 newnorm[2];
    newnorm[0] = norm * cangle - norfpc0 * sangle;
    newnorm[1] = -norm * cangle - norfpc0 * sangle;
    newnorm[0].normalize();
    newnorm[1].normalize();

    st0.v[0] = c0 + r0 * newnorm[0];
    st0.v[1] = c1 + r1 * newnorm[0];
    st0.v[2] = c2 + r2 * newnorm[0];
    // st0.update_normal();
    st0.update_normal(newnorm[0]);

    st1.v[0] = c0 + r0 * newnorm[1];
    st1.v[1] = c1 + r1 * newnorm[1];
    st1.v[2] = c2 + r2 * newnorm[1];
    // st1.update_normal(true);
    st1.update_normal(newnorm[1]);
  }

  return type;
}

}  // namespace matfp