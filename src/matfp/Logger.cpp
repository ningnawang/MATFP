// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include <matfp/Logger.h>
#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <memory>
#include <vector>

namespace matfp {

std::shared_ptr<spdlog::logger> Logger::logger_;

// Some code was copied over from <spdlog/async.h>
void Logger::init(std::string log_name, bool use_cout,
                  const std::string &filename, bool truncate) {
  std::vector<spdlog::sink_ptr> sinks;
  if (use_cout) {
    sinks.emplace_back(std::make_shared<spdlog::sinks::stdout_color_sink_mt>());
  }
  if (!filename.empty()) {
    sinks.emplace_back(std::make_shared<spdlog::sinks::basic_file_sink_mt>(
        filename, truncate));
  }

  spdlog::init_thread_pool(8192, 1);

  logger_ = std::make_shared<spdlog::async_logger>(
      log_name, sinks.begin(), sinks.end(), spdlog::thread_pool(),
      spdlog::async_overflow_policy::block);

  spdlog::drop(log_name);
  spdlog::register_logger(logger_);
}

}  // namespace matfp
