// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>
#include <matfp/Logger.h>
#include <spdlog/details/registry.h>
#include <spdlog/details/thread_pool.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <iostream>
#include <memory>
#include <mutex>

namespace matfp {

std::shared_ptr<spdlog::async_logger> Logger::logger_;

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

  auto &registry_inst = spdlog::details::registry::instance();

  // create global thread pool if not already exists..
  std::lock_guard<std::recursive_mutex> tp_lock(registry_inst.tp_mutex());
  auto tp = registry_inst.get_tp();
  if (tp == nullptr) {
    tp = std::make_shared<spdlog::details::thread_pool>(
        spdlog::details::default_async_q_size, 1);
    registry_inst.set_tp(tp);
  }

  logger_ = std::make_shared<spdlog::async_logger>(
      log_name, sinks.begin(), sinks.end(), std::move(tp),
      spdlog::async_overflow_policy::block);
  registry_inst.register_and_init(logger_);
}

}  // namespace matfp
