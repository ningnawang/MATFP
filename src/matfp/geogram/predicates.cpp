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
#include <geogram/basic/assert.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/common.h>
#include <geogram/basic/logger.h>
#include <geogram/basic/matrix.h>
#include <geogram/numerics/multi_precision.h>
#include <geogram/numerics/predicates.h>
#include <matfp/Logger.h>
#include <matfp/Types/CommonTypes.h>
#include <matfp/geogram/predicates.h>

#include <algorithm>

#define FPG_UNCERTAIN_VALUE 0
#include <matfp/geogram/predicates/powerside1.h>
#include <matfp/geogram/predicates/powerside2.h>
#include <matfp/geogram/predicates/powerside3.h>
#include <matfp/geogram/predicates/powerside4.h>

// Reference: https://hal.inria.fr/hal-01225202/document
namespace matfp {
using namespace GEO;
using namespace std;

// ninwang 2021-03-08:
// not using GEO::PCK::SOS_ADDRESS here
// this might becaues our [double* p0] has variable length?
// it store weight at 4th position as well
// and that made Symbolic perturbation really sad :(
GEO::PCK::SOSMode SOS_mode_ = GEO::PCK::SOS_LEXICO;

/**
 * \brief Compares two 3D points with respect to the lexicographic
 *  order.
 * \param[in] x , y pointers to the coordinates of the two 3D points.
 * \retval true if x is strictly before y in the lexicographic order.
 * \retval false otherwise.
 */
bool lexico_compare_3d(const double* x, const double* y) {
  if (x[0] < y[0]) {
    return true;
  }
  if (x[0] > y[0]) {
    return false;
  }
  if (x[1] < y[1]) {
    return true;
  }
  if (x[1] > y[1]) {
    return false;
  }
  return x[2] < y[2];
}

/**
 * \brief Sorts an array of pointers to points.
 * \details set_SOS_mode() alters the behavior of this function.
 *  If set to PCK::SOS_ADDRESS, then just the addresses of the points
 *  are sorted. If set to PCK::SOS_LEXICO, then the points are sorted
 *  in function of the lexicographic order of their coordinates.
 * \param[in] begin a pointer to the first point.
 * \param[in] end one position past the pointer to the last point.
 * \param[in] dim the dimension of the points.
 */
void SOS_sort(const double** begin, const double** end, index_t dim) {
  if (SOS_mode_ == GEO::PCK::SOS_ADDRESS) {
    std::sort(begin, end);
  } else {
    // if(dim == 3) {
    std::sort(begin, end, lexico_compare_3d);
    // } else {
    //     std::sort(begin, end, LexicoCompare(dim));
    // }
  }
}

// ================= side1 =========================================

/**
 * \brief Exact implementation of the side1() predicate using low-level
 *  exact arithmetics API (expansion class).
 */
GEO::Sign side1_exact_SOS(const double* p0, const double w0, const double* p1,
                          const double w1, const double* q0,
                          coord_index_t dim) {
  expansion& l = expansion_sq_dist(p0, p1, dim);
  expansion& R = expansion_diff(w0, w1);
  expansion& L = expansion_sum(l, R);

  expansion& a = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
  expansion& r = expansion_diff(L, a);
  Sign r_sign = r.sign();

  // logger().debug("r_sign in side1_exact_SOS: {} ", r_sign);

  // Symbolic perturbation, Simulation of Simplicity
  if (r_sign == ZERO) {
    // logger().debug("p0: ({},{},{},{}), p1 ({},{},{},{}), q0 ({},{},{})",
    //     p0[0], p0[1], p0[2], p0[3],
    //     p1[0], p1[1], p1[2], p1[3],
    //     q0[0], q0[1], q0[2]
    // );

    // int is_positive = p0 < p1 ?  1 : -1;
    // logger().debug("p0 < p1 ? {}", is_positive);
    // return (p0 < p1) ? POSITIVE : NEGATIVE;

    /////////
    // why p0 < p1 not giving consistant result
    // for the exact same input???
    // ninwang: i guess the length of double* also affects the results
    for (index_t i = 0; i < dim; ++i) {
      if (p0[i] < p1[i]) return POSITIVE;
      if (p0[i] > p1[i]) return NEGATIVE;
    }
    geo_assert_not_reached;
  }
  return r_sign;
}
// ================= side2 =========================================

/**
 * \brief Exact implementation of the side2() predicate using low-level
 *  exact arithmetics API (expansion class).
 */
GEO::Sign side2_exact_SOS(const double* p0, const double w0, const double* p1,
                          const double w1, const double* p2, const double w2,
                          const double* q0, const double* q1,
                          GEO::coord_index_t dim) {
  const expansion& l1 = expansion_sq_dist(p1, p0, dim);
  const expansion& l2 = expansion_sq_dist(p2, p0, dim);

  const expansion& R1 = expansion_diff(w0, w1);
  const expansion& R2 = expansion_diff(w0, w2);

  const expansion& L1 = expansion_sum(l1, R1);
  const expansion& L2 = expansion_sum(l2, R2);

  const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
  const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
  const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
  const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);

  const expansion& Delta = expansion_diff(a11, a10);

  Sign Delta_sign = Delta.sign();
  // Should not occur with symbolic
  // perturbation done at previous steps.
  geo_assert(Delta_sign != ZERO);

  //       [ Lambda0 ]   [ -1 ]        [  a11 ]
  // Delta [         ] = [    ] * L1 + [      ]
  //       [ Lambda1 ]   [  1 ]        [ -a10 ]

  const expansion& DeltaLambda0 = expansion_diff(a11, L1);
  const expansion& DeltaLambda1 = expansion_diff(L1, a10);

  // r = Delta*L2 - ( a20*DeltaLambda0 + a21*DeltaLambda1 )

  const expansion& r0 = expansion_product(Delta, L2);
  const expansion& r1 = expansion_product(a20, DeltaLambda0).negate();
  const expansion& r2 = expansion_product(a21, DeltaLambda1).negate();
  const expansion& r = expansion_sum3(r0, r1, r2);

  Sign r_sign = r.sign();

  // Simulation of Simplicity (symbolic perturbation)
  if (r_sign == ZERO) {
    const double* p_sort[3];
    p_sort[0] = p0;
    p_sort[1] = p1;
    p_sort[2] = p2;

    SOS_sort(p_sort, p_sort + 3, dim);

    for (index_t i = 0; i < 3; ++i) {
      if (p_sort[i] == p0) {
        const expansion& z1 = expansion_diff(Delta, a21);
        const expansion& z = expansion_sum(z1, a20);
        Sign z_sign = z.sign();
        // len_side2_SOS = std::max(len_side2_SOS, z.length());
        if (z_sign != ZERO) {
          return Sign(Delta_sign * z_sign);
        }
      }
      if (p_sort[i] == p1) {
        const expansion& z = expansion_diff(a21, a20);
        Sign z_sign = z.sign();
        // len_side2_SOS = std::max(len_side2_SOS, z.length());
        if (z_sign != ZERO) {
          return Sign(Delta_sign * z_sign);
        }
      }
      if (p_sort[i] == p2) {
        return NEGATIVE;
      }
    }
    geo_assert_not_reached;
  }

  return Sign(Delta_sign * r_sign);
}

// ================= side3 =========================================

/**
 * \brief Exact implementation of the side3() predicate using low-level
 *  exact arithmetics API (expansion class).
 */
Sign side3_exact_SOS(const double* p0, const double w0, const double* p1,
                     const double w1, const double* p2, const double w2,
                     const double* p3, const double w3, const double* q0,
                     const double* q1, const double* q2, coord_index_t dim) {
  const expansion& l1 = expansion_sq_dist(p1, p0, dim);
  const expansion& l2 = expansion_sq_dist(p2, p0, dim);
  const expansion& l3 = expansion_sq_dist(p3, p0, dim);

  const expansion& R1 = expansion_diff(w0, w1);
  const expansion& R2 = expansion_diff(w0, w2);
  const expansion& R3 = expansion_diff(w0, w3);

  const expansion& L1 = expansion_sum(l1, R1);
  const expansion& L2 = expansion_sum(l2, R2);
  const expansion& L3 = expansion_sum(l3, R3);

  const expansion& a10 = expansion_dot_at(p1, q0, p0, dim).scale_fast(2.0);
  const expansion& a11 = expansion_dot_at(p1, q1, p0, dim).scale_fast(2.0);
  const expansion& a12 = expansion_dot_at(p1, q2, p0, dim).scale_fast(2.0);
  const expansion& a20 = expansion_dot_at(p2, q0, p0, dim).scale_fast(2.0);
  const expansion& a21 = expansion_dot_at(p2, q1, p0, dim).scale_fast(2.0);
  const expansion& a22 = expansion_dot_at(p2, q2, p0, dim).scale_fast(2.0);

  const expansion& a30 = expansion_dot_at(p3, q0, p0, dim).scale_fast(2.0);
  const expansion& a31 = expansion_dot_at(p3, q1, p0, dim).scale_fast(2.0);
  const expansion& a32 = expansion_dot_at(p3, q2, p0, dim).scale_fast(2.0);

  // [ b00 b01 b02 ]           [  1   1   1  ]-1
  // [ b10 b11 b12 ] = Delta * [ a10 a11 a12 ]
  // [ b20 b21 b22 ]           [ a20 a21 a22 ]

  const expansion& b00 = expansion_det2x2(a11, a12, a21, a22);
  const expansion& b01 = expansion_diff(a21, a22);
  const expansion& b02 = expansion_diff(a12, a11);
  const expansion& b10 = expansion_det2x2(a12, a10, a22, a20);
  const expansion& b11 = expansion_diff(a22, a20);
  const expansion& b12 = expansion_diff(a10, a12);
  const expansion& b20 = expansion_det2x2(a10, a11, a20, a21);
  const expansion& b21 = expansion_diff(a20, a21);
  const expansion& b22 = expansion_diff(a11, a10);

  const expansion& Delta = expansion_sum3(b00, b10, b20);
  Sign Delta_sign = Delta.sign();
  // Should not occur with symbolic
  // perturbation done at previous steps.
  geo_assert(Delta_sign != ZERO);

  //       [ Lambda0 ]   [ b01 b02 ]   [ L1 ]   [ b00 ]
  // Delta [ Lambda1 ] = [ b11 b12 ] * [    ] + [ b10 ]
  //       [ Lambda2 ]   [ b21 b22 ]   [ L2 ]   [ b20 ]

  const expansion& b01_l1 = expansion_product(b01, L1);
  const expansion& b02_l2 = expansion_product(b02, L2);
  const expansion& DeltaLambda0 = expansion_sum3(b01_l1, b02_l2, b00);

  const expansion& b11_l1 = expansion_product(b11, L1);
  const expansion& b12_l2 = expansion_product(b12, L2);
  const expansion& DeltaLambda1 = expansion_sum3(b11_l1, b12_l2, b10);

  const expansion& b21_l1 = expansion_product(b21, L1);
  const expansion& b22_l2 = expansion_product(b22, L2);
  const expansion& DeltaLambda2 = expansion_sum3(b21_l1, b22_l2, b20);

  // r = Delta*L3-(a30*DeltaLambda0+a31*DeltaLambda1+a32*DeltaLambda2)

  const expansion& r0 = expansion_product(Delta, L3);
  const expansion& r1 = expansion_product(a30, DeltaLambda0).negate();
  const expansion& r2 = expansion_product(a31, DeltaLambda1).negate();
  const expansion& r3 = expansion_product(a32, DeltaLambda2).negate();
  const expansion& r = expansion_sum4(r0, r1, r2, r3);
  Sign r_sign = r.sign();

  // Simulation of Simplicity (symbolic perturbation)
  if (r_sign == ZERO) {
    const double* p_sort[4];
    p_sort[0] = p0;
    p_sort[1] = p1;
    p_sort[2] = p2;
    p_sort[3] = p3;
    SOS_sort(p_sort, p_sort + 4, dim);
    for (index_t i = 0; i < 4; ++i) {
      if (p_sort[i] == p0) {
        const expansion& z1_0 = expansion_sum(b01, b02);
        const expansion& z1 = expansion_product(a30, z1_0).negate();
        const expansion& z2_0 = expansion_sum(b11, b12);
        const expansion& z2 = expansion_product(a31, z2_0).negate();
        const expansion& z3_0 = expansion_sum(b21, b22);
        const expansion& z3 = expansion_product(a32, z3_0).negate();
        const expansion& z = expansion_sum4(Delta, z1, z2, z3);
        Sign z_sign = z.sign();
        if (z_sign != ZERO) {
          return Sign(Delta_sign * z_sign);
        }
      } else if (p_sort[i] == p1) {
        const expansion& z1 = expansion_product(a30, b01);
        const expansion& z2 = expansion_product(a31, b11);
        const expansion& z3 = expansion_product(a32, b21);
        const expansion& z = expansion_sum3(z1, z2, z3);
        Sign z_sign = z.sign();
        if (z_sign != ZERO) {
          return Sign(Delta_sign * z_sign);
        }
      } else if (p_sort[i] == p2) {
        const expansion& z1 = expansion_product(a30, b02);
        const expansion& z2 = expansion_product(a31, b12);
        const expansion& z3 = expansion_product(a32, b22);
        const expansion& z = expansion_sum3(z1, z2, z3);
        Sign z_sign = z.sign();
        if (z_sign != ZERO) {
          return Sign(Delta_sign * z_sign);
        }
      } else if (p_sort[i] == p3) {
        return NEGATIVE;
      }
    }
    geo_assert_not_reached;
  }
  return Sign(Delta_sign * r_sign);
}

// ================= side4 =========================================
/**
 * \brief Exact implementation of the side4_3d_SOS() predicate
 *  using low-level exact arithmetics API (expansion class).
 * \param[in] sos if true, applies symbolic perturbation when
 *  result is zero, else returns zero
 */
Sign power_side4_3d_exact_SOS(const double* p0, const double w0,
                              const double* p1, const double w1,
                              const double* p2, const double w2,
                              const double* p3, const double w3,
                              const double* p4, const double w4,
                              bool sos = true) {
  const expansion& l1 = expansion_sq_dist(p1, p0, 3);
  const expansion& l2 = expansion_sq_dist(p2, p0, 3);
  const expansion& l3 = expansion_sq_dist(p3, p0, 3);
  const expansion& l4 = expansion_sq_dist(p4, p0, 3);

  const expansion& R1 = expansion_diff(w0, w1);
  const expansion& R2 = expansion_diff(w0, w2);
  const expansion& R3 = expansion_diff(w0, w3);
  const expansion& R4 = expansion_diff(w0, w4);

  const expansion& a11 = expansion_diff(p1[0], p0[0]);
  const expansion& a12 = expansion_diff(p1[1], p0[1]);
  const expansion& a13 = expansion_diff(p1[2], p0[2]);
  const expansion& a14 = expansion_sum(l1, R1).negate();

  const expansion& a21 = expansion_diff(p2[0], p0[0]);
  const expansion& a22 = expansion_diff(p2[1], p0[1]);
  const expansion& a23 = expansion_diff(p2[2], p0[2]);
  const expansion& a24 = expansion_sum(l2, R2).negate();

  const expansion& a31 = expansion_diff(p3[0], p0[0]);
  const expansion& a32 = expansion_diff(p3[1], p0[1]);
  const expansion& a33 = expansion_diff(p3[2], p0[2]);
  const expansion& a34 = expansion_sum(l3, R3).negate();

  const expansion& a41 = expansion_diff(p4[0], p0[0]);
  const expansion& a42 = expansion_diff(p4[1], p0[1]);
  const expansion& a43 = expansion_diff(p4[2], p0[2]);
  const expansion& a44 = expansion_sum(l4, R4).negate();

  // This commented-out version does not reuse
  // the 2x2 minors.
  /*
          const expansion& Delta1 = expansion_det3x3(
              a21, a22, a23,
              a31, a32, a33,
              a41, a42, a43
          );
          const expansion& Delta2 = expansion_det3x3(
              a11, a12, a13,
              a31, a32, a33,
              a41, a42, a43
          );
          const expansion& Delta3 = expansion_det3x3(
              a11, a12, a13,
              a21, a22, a23,
              a41, a42, a43
          );
          const expansion& Delta4 = expansion_det3x3(
              a11, a12, a13,
              a21, a22, a23,
              a31, a32, a33
          );
  */

  // Optimized version that reuses the 2x2 minors

  const expansion& m12 = expansion_det2x2(a12, a13, a22, a23);
  const expansion& m13 = expansion_det2x2(a12, a13, a32, a33);
  const expansion& m14 = expansion_det2x2(a12, a13, a42, a43);
  const expansion& m23 = expansion_det2x2(a22, a23, a32, a33);
  const expansion& m24 = expansion_det2x2(a22, a23, a42, a43);
  const expansion& m34 = expansion_det2x2(a32, a33, a42, a43);

  const expansion& z11 = expansion_product(a21, m34);
  const expansion& z12 = expansion_product(a31, m24).negate();
  const expansion& z13 = expansion_product(a41, m23);
  const expansion& Delta1 = expansion_sum3(z11, z12, z13);

  const expansion& z21 = expansion_product(a11, m34);
  const expansion& z22 = expansion_product(a31, m14).negate();
  const expansion& z23 = expansion_product(a41, m13);
  const expansion& Delta2 = expansion_sum3(z21, z22, z23);

  const expansion& z31 = expansion_product(a11, m24);
  const expansion& z32 = expansion_product(a21, m14).negate();
  const expansion& z33 = expansion_product(a41, m12);
  const expansion& Delta3 = expansion_sum3(z31, z32, z33);

  const expansion& z41 = expansion_product(a11, m23);
  const expansion& z42 = expansion_product(a21, m13).negate();
  const expansion& z43 = expansion_product(a31, m12);
  const expansion& Delta4 = expansion_sum3(z41, z42, z43);

  Sign Delta4_sign = Delta4.sign();
  geo_assert(Delta4_sign != ZERO);

  const expansion& r_1 = expansion_product(Delta1, a14);
  const expansion& r_2 = expansion_product(Delta2, a24).negate();
  const expansion& r_3 = expansion_product(Delta3, a34);
  const expansion& r_4 = expansion_product(Delta4, a44).negate();
  const expansion& r = expansion_sum4(r_1, r_2, r_3, r_4);
  Sign r_sign = r.sign();

  // Simulation of Simplicity (symbolic perturbation)
  if (sos && r_sign == ZERO) {
    const double* p_sort[5];
    p_sort[0] = p0;
    p_sort[1] = p1;
    p_sort[2] = p2;
    p_sort[3] = p3;
    p_sort[4] = p4;
    SOS_sort(p_sort, p_sort + 5, 3);
    for (index_t i = 0; i < 5; ++i) {
      if (p_sort[i] == p0) {
        const expansion& z1 = expansion_diff(Delta2, Delta1);
        const expansion& z2 = expansion_diff(Delta4, Delta3);
        const expansion& z = expansion_sum(z1, z2);
        Sign z_sign = z.sign();
        if (z_sign != ZERO) {
          return Sign(Delta4_sign * z_sign);
        }
      } else if (p_sort[i] == p1) {
        Sign Delta1_sign = Delta1.sign();
        if (Delta1_sign != ZERO) {
          return Sign(Delta4_sign * Delta1_sign);
        }
      } else if (p_sort[i] == p2) {
        Sign Delta2_sign = Delta2.sign();
        if (Delta2_sign != ZERO) {
          return Sign(-Delta4_sign * Delta2_sign);
        }
      } else if (p_sort[i] == p3) {
        Sign Delta3_sign = Delta3.sign();
        if (Delta3_sign != ZERO) {
          return Sign(Delta4_sign * Delta3_sign);
        }
      } else if (p_sort[i] == p4) {
        return NEGATIVE;
      }
    }
  }
  return Sign(Delta4_sign * r_sign);
}

inline void print_sign(Sign s) {
  switch (s) {
    case Sign::ZERO:
      std::cout << "side3_3d_SOS sign is 0\n";
      break;
    case Sign::POSITIVE:
      std::cout << "side3_3d_SOS sign is 1\n";
      break;
    case Sign::NEGATIVE:
      std::cout << "side3_3d_SOS sign is -1\n";
      break;
    default:
      break;
  }
}

}  // namespace matfp

namespace matfp {

GEO::Sign PCK::power_side1_SOS(const double* p0, const double w0,
                               const double* p1, const double w1,
                               const double* q0) {
  GEO::Sign result;
  result = GEO::Sign(side1_power_3d_filter(p0, w0, p1, w1, q0));
  // logger().debug("sign of side1_power_3d_filter: {}", result);
  if (result == GEO::ZERO) {
    result = GEO::Sign(side1_exact_SOS(p0, w0, p1, w1, q0, 3));
  }
  return result;
}

GEO::Sign PCK::power_side2_SOS(const double* p0, const double w0,
                               const double* p1, const double w1,
                               const double* p2, const double w2,
                               const double* q0, const double* q1) {
  GEO::Sign result;
  result = GEO::Sign(side2_3d_filter(p0, w0, p1, w1, p2, w2, q0, q1));
  if (result == GEO::ZERO) {
    result = GEO::Sign(side2_exact_SOS(p0, w0, p1, w1, p2, w2, q0, q1, 3));
  }
  return result;
}

GEO::Sign PCK::power_side3_SOS(const double* p0, const double w0,
                               const double* p1, const double w1,
                               const double* p2, const double w2,
                               const double* p3, const double w3,
                               const double* q0, const double* q1,
                               const double* q2) {
  GEO::Sign result;
  result =
      GEO::Sign(side3_3d_filter(p0, w0, p1, w1, p2, w2, p3, w3, q0, q1, q2));
  if (result == GEO::ZERO) {
    result = GEO::Sign(
        side3_exact_SOS(p0, w0, p1, w1, p2, w2, p3, w3, q0, q1, q2, 3));
  }
  return result;
}

GEO::Sign PCK::power_side4_3d_SOS(const double* p0, const double w0,
                                  const double* p1, const double w1,
                                  const double* p2, const double w2,
                                  const double* p3, const double w3,
                                  const double* p4, const double w4) {
  GEO::Sign result;
  result = GEO::Sign(side4_3d_filter(p0, w0, p1, w1, p2, w2, p3, w3, p4, w4));
  if (result == ZERO) {
    result = power_side4_3d_exact_SOS(p0, w0, p1, w1, p2, w2, p3, w3, p4, w4);
  }
  return result;
}

// GEO::Sign PCK::power_side4_SOS(
//     const double* p0, const double w0,
//     const double* p1, const double w1,
//     const double* p2, const double w2,
//     const double* p3, const double w3,
//     const double* p4, const double w4,
//     const double* q0, const double* q1,
//     const double* q2, const double* q3
// ) {
//     GEO::Sign result;
//     result = GEO::Sign(side4_3d_filter(p0, w0, p1, w1, p2, w2, p3, w3, p4,
//     w4, q0, q1, q2, q3)); if(result == ZERO) {
//         result = power_side4_3d_exact_SOS(p0, w0, p1, w1, p2, w2, p3, w3, p4,
//         w4);
//     }
//     return result;
// }
}  // namespace matfp
