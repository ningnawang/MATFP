#pragma once

#include <fmt/ranges.h>
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>
#include <matfp/Exception.h>
#include <spdlog/async.h>
#include <spdlog/spdlog.h>

namespace matfp {

struct Logger {
  static std::shared_ptr<spdlog::async_logger> logger_;

  // By default, write to stdout, but don't write to any file
  static void init(std::string log_name = "", bool use_cout = true,
                   const std::string &filename = "", bool truncate = true);
};

// Retrieve current logger, or create one if not available
inline spdlog::async_logger &logger() {
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
