// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#pragma once

// Must come first
#include <matfp/LoggerFormatter.h>

// These can come after
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>
#include <matfp/Exception.h>
#include <spdlog/spdlog.h>

#include <memory>

namespace matfp {

struct Logger {
  static std::shared_ptr<spdlog::logger> logger_;

  static void init(std::string log_name = "", bool use_cout = true,
                   const spdlog::filename_t& filename = spdlog::filename_t(),
                   bool truncate = true);
};

inline spdlog::logger& logger() {
  if (!Logger::logger_) {
    Logger::init();
  }
  return *Logger::logger_;
}

template <typename T>
[[noreturn]] void log_and_throw(T x) {
  logger().error(x);
  throw MATFPError(x);
}

}  // namespace matfp
