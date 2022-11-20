// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>

#include <Eigen/Dense>
#include <limits>
#include <queue>

namespace matfp {

enum SphereType {
  T_UNK = -1,
  T_1_N = 1,  // mesh feature, sharp edges, and corners (abuse)
  T_2 = 2,
  T_N = 3,     // including T_3, T_4 ...
  T_2_c = 11,  // added for one concave line using sphere shrinking
  T_N_c = 12   // added through internal feature preservation
};

enum DeletionType {
  D_UNK = -1,    // unknown
  D_NOT = 0,     // not deleted
  D_FP_EXT = 1,  // deleted for FP external deletion
  D_CLEAN = 3,   // deleted for cleaning
  D_DUP = 4      // duplications
};

enum AdditionType {
  A_NOT = 0,     // not sampled
  A_FP_EXT = 1,  // sampled for FP external corner
  A_FP_INT = 2   // sampled for FP internal addition
};

// degree between two normals
// if less than 10, then we think they are on the same plane
// 5 degree is too small, for two vertices on a curved face
// their degree might be > 5
const double esp_degree_5 = 5;
const double esp_degree_10 = 10;
const double esp_degree_20 = 20;
const double esp_degree_30 = 30;

// for Feature Distance
const int FEATURE_FD_RHO = -1;
const size_t FEATURE_FD = size_t(-1);

// constants
const double PI = M_PI;
const double HALF_PI = M_PI / 180.;
const double STD_NAN = std::numeric_limits<double>::quiet_NaN();

#ifdef MATFP_USE_FLOAT
typedef float Scalar;
#define SCALAR_ZERO 1e-6
#define SCALAR_ZERO_2 1e-12
#define SCALAR_ZERO_3 1e-18
#else
typedef double Scalar;
#define SCALAR_ZERO_N3 1e-3
#define SCALAR_ZERO_N4 1e-4
#define SCALAR_ZERO_N5 1e-5
#define SCALAR_ZERO_N6 1e-6
#define SCALAR_ZERO 1e-8
#define SCALAR_ZERO_2 1e-16
#define SCALAR_ZERO_3 1e-24
#define SCALAR_FEATURE_SQ_RADIUS 1e-24
#endif

const double extrude_threshold = 0.05;
const double extrude_threshold_2 = 0.01;
const double INIT_RADIUS = 10;  // for sphere shriking algo
const double SQ_DIST_TO_CC =
    SCALAR_ZERO_N5;  // if sphere too close to tan_cc_lines, then delete

///////
// reversed priority queue, smaller first
// C++11 added alias declarations, which are generalization of typedef
template <typename T>
using MIN_PQ = std::priority_queue<T, std::vector<T>, std::greater<T>>;

////////////
// Eigen
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> VectorXs;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<Scalar, 2, 2> Matrix2;
typedef Eigen::Matrix<Scalar, 3, 3> Matrix3;
typedef Eigen::Matrix<Scalar, 4, 4> Matrix4;
typedef Eigen::Matrix<Scalar, 4, 1> Vector4;
typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
typedef Eigen::Matrix<Scalar, 2, 1> Vector2;
typedef Eigen::Matrix<int, 4, 1> Vector4i;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<int, 2, 1> Vector2i;

//////////
// CGAL
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kf;
typedef Kf::FT Weight;
typedef Kf::Point_3 Point;
typedef Kf::Weighted_point_3 Weighted_point;

}  // namespace matfp
