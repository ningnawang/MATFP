// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

#include <stdexcept>
#include <string>

namespace matfp {

class MATFPError : public std::runtime_error {
 public:
  explicit MATFPError(const std::string& what_arg)
      : std::runtime_error(what_arg) {}
};

}  // namespace matfp
