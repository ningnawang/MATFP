// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <geogram/mesh/mesh.h>

#include <Eigen/Dense>
#include <unordered_set>

namespace pre_matfp {

class Shape3D {
 public:
  GEO::Mesh sf_mesh;
};

#define SCALAR_ZERO 1e-8
const double PI = M_PI;

////////////
// Eigen
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXs;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXs;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> VectorXi;
typedef Eigen::Matrix<double, 3, 3> Matrix3;
typedef Eigen::Matrix<double, 3, 1> Vector3;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<int, 2, 1> Vector2i;

inline Vector3 to_eigen(const GEO::vec3& v) {
  return Vector3(v[0], v[1], v[2]);
}

template <typename T>
void vector_unique(std::vector<T>& v) {
  std::sort(v.begin(), v.end());
  v.erase(std::unique(v.begin(), v.end()), v.end());
}

inline Vector3 get_normal(const Vector3& a, const Vector3& b,
                          const Vector3& c) {
  // facet abc is clockwise/counter-clockwise
  // (a-c) x (b-c) is clockwise/counter-clockwise
  return ((a - c).cross(b - c)).normalized();
}

inline Vector3 get_triangle_centroid(const Vector3& v1, const Vector3& v2,
                                     const Vector3& v3) {
  return (v1 + v2 + v3) / 3.;
};

template <typename T>
void set_intersection(const std::unordered_set<T>& s1,
                      const std::unordered_set<T>& s2, std::vector<T>& v) {
  if (s2.size() < s1.size()) {
    set_intersection(s2, s1, v);
    return;
  }
  v.clear();
  v.reserve(std::min(s1.size(), s2.size()));
  for (int x : s1) {
    if (s2.count(x)) {
      v.push_back(x);
    }
  }
}

inline int mod3(int j) { return j % 3; }

}  // namespace pre_matfp
