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
#include <geogram/basic/numeric.h>

namespace matfp {

// Reference: https://hal.inria.fr/hal-01225202/document
// DIM = 3 ONLY, for power diagram
namespace PCK {
/**
 * \brief Computes the side of a point (given directly)
 *  relative to a power bisector.
 * \details Computes the side of \f$ q0 \f$ relative to
 * \f$ \Pi(p0,p1) \f$.
 * Symbolic perturbation is applied whenever equality holds.
 * \param[in] p0 , p1 extremities of the power bisector
 * \param[in] q0 point to be tested
 * \param[in] DIM number of coordinates of the point
 * \retval POSITIVE if d(p0,q0) < d(p1,q0)
 * \retval NEGATIVE if d(p0,q0) > d(p1,q1)
 * \retval perturb() if f(p0,q0) = d(p1,q1),
 *  where \c perturb() denotes a globally
 *  consistent perturbation, that returns either POSITIVE or NEGATIVE
 * \note Only some specific dimensions are implemented (3,4,6 and 7)
 */
GEO::Sign power_side1_SOS(const double* p0, const double* p1, const double* q0);

GEO::Sign power_side1_SOS(const double* p0, const double w0, const double* p1,
                          const double w1, const double* q0);

/**
 * \brief Computes the side of a point (given as the intersection
 *  between a segment and a power bisector) relative to another power bisector.
 * \details Computes the side of \f$ q = \Pi(p0,p1) \cap [q0,q1] \f$
 * relative to \f$ \Pi(p0,p2) \f$.
 * Symbolic perturbation is applied whenever equality holds.
 * \param[in] p0 first extremity of the power bisectors
 * \param[in] p1 second extremity of the first power bisector
 *  (that defines the intersection q)
 * \param[in] p2 second extremity of the second power bisector
 *  (against which orientation is tested)
 * \param[in] q0 , q1 extremities of the segment
 *  (that defines the intersection q)
 * \retval POSITIVE if d(p0,q) < d(p2,q)
 * \retval NEGATIVE if d(p0,q) > d(p2,q)
 * \retval perturb() if d(p0,q) = d(p2,q),
 *  where \c perturb() denotes a globally
 *  consistent perturbation, that returns either POSITIVE or NEGATIVE
 * \note Only some specific dimensions are implemented (3,4,6 and 7)
 */
GEO::Sign power_side2_SOS(const double* p0, const double* p1, const double* p2,
                          const double* q0, const double* q1);
GEO::Sign power_side2_SOS(const double* p0, const double w0, const double* p1,
                          const double w1, const double* p2, const double w2,
                          const double* q0, const double* q1);

/**
 * \brief Computes the side of a point (given as the intersection
 *  between a facet and two power bisectors) relative to another power bisector.
 * \details Computes the side of
 *  \f$ q = \Pi(p0,p1) \cap Pi(p0,p2) \cap \Delta[q0,q1,q2] \f$
 * relative to \f$ \Pi(p0,p3) \f$.
 * Symbolic perturbation is applied whenever equality holds.
 * \param[in] p0 first extremity of the power bisectors
 * \param[in] p1 second extremity of the first power bisector
 *  (that defines the intersection q)
 * \param[in] p2 second extremity of the second power bisector
 *  (that defines the intersection q)
 * \param[in] p3 second extremity of the third power bisector
 *  (against which orientation is tested)
 * \param[in] q0 , q1 , q2 vertices of the triangle
 *  (that defines the intersection q)
 * \retval POSITIVE if d(p0,q) < d(p3,q)
 * \retval NEGATIVE if d(p0,q) > d(p3,q)
 * \retval perturb() if d(p0,q) = d(p3,q),
 *  where \c perturb() denotes a globally
 *  consistent perturbation, that returns either POSITIVE or NEGATIVE
 * \note Only some specific dimensions are implemented (3,4,6 and 7)
 */
GEO::Sign power_side3_SOS(const double* p0, const double* p1, const double* p2,
                          const double* p3, const double* q0, const double* q1,
                          const double* q2);
GEO::Sign power_side3_SOS(const double* p0, const double w0, const double* p1,
                          const double w1, const double* p2, const double w2,
                          const double* p3, const double w3, const double* q0,
                          const double* q1, const double* q2);

/**
 * \brief Computes the side of a point (given as the intersection
 *   between three bisectors) relative to another bisector.
 * \details Computes the side of
 *  \f$ q = \Pi(p0,p1) \cap \Pi(p0,p2) \cap \Pi(p0,p3) \f$
 * relative to \f$ Pi(p0,p4) \f$.
 * Symbolic perturbation is applied whenever equality holds.
 * side4_3d() is a special case of side4(), where the ambient and
 * intrinsic dimensions coincide (therefore no embedding tetrahedron
 * is needed).
 * \param[in] p0 first extremity of the bisectors
 * \param[in] p1 second extremity of the first bisector
 *  (that defines the intersection q)
 * \param[in] p2 second extremity of the second bisector
 *  (that defines the intersection q)
 * \param[in] p3 second extremity of the third bisector
 *  (that defines the intersection q)
 * \param[in] p4 second extremity of the fourth bisector
 *  (against which orientation is tested)
 * \retval POSITIVE if d(p0,q) < d(p4,q)
 * \retval NEGATIVE if d(p0,q) > d(p4,q)
 * \retval perturb() if d(p0,q) = d(p4,q),
 *  where \c perturb() denotes a globally
 *  consistent perturbation, that returns either POSITIVE or NEGATIVE
 */
GEO::Sign power_side4_3d_SOS(const double* p0, const double w0,
                             const double* p1, const double w1,
                             const double* p2, const double w2,
                             const double* p3, const double w3,
                             const double* p4, const double w4);

// /**
//  * \brief Computes the side of a point (given as the intersection
//  *   between a tetrahedron and three bisectors) relative to
//  *  another bisector.
//  * \details Computes the side of
//  *  \f$ q = \Pi(p0,p1) \cap Pi(p0,p2) \cap Pi(p0,p3)
//  * \cap \Delta[q0,q1,q2,q3] \f$ relative to \f$ \Pi(p0,p4) \f$.
//  * Symbolic perturbation is applied whenever equality holds.
//  * \param[in] p0 first extremity of the bisectors
//  * \param[in] p1 second extremity of the first bisector
//  *  (that defines the intersection q)
//  * \param[in] p2 second extremity of the second bisector
//  *  (that defines the intersection q)
//  * \param[in] p3 second extremity of the third bisector
//  *  (that defines the intersection q)
//  * \param[in] p4 second extremity of the fourth bisector
//  *  (against which orientation is tested)
//  * \param[in] q0 , q1 , q2 , q3 vertices of the tetrahedron
//  *  (that defines the intersection q)
//  *  (that defines the intersection q)
//  * \retval POSITIVE if d(p0,q) < d(p4,q)
//  * \retval NEGATIVE if d(p0,q) > d(p4,q)
//  * \retval perturb() if d(p0,q) = d(p4,q),
//  *  where \c perturb() denotes a globally
//  *  consistent perturbation, that returns either POSITIVE or NEGATIVE
//  * \note Only some specific dimensions are implemented (3,4,6 and 7)
//  */
// GEO::Sign power_side4_SOS(
//     const double* p0, const double w0,
//     const double* p1, const double w1,
//     const double* p2, const double w2,
//     const double* p3, const double w3,
//     const double* p4, const double w4,
//     const double* q0 = nullptr, const double* q1 = nullptr,
//     const double* q2 = nullptr, const double* q3 = nullptr
// );

}  // namespace PCK

}  // namespace  matfp
