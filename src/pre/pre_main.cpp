// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <matfp/Args.h>
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>
#include <matfp/Logger.h>

#include <CLI/CLI.hpp>

#include "pre_guiwindow.h"

using namespace matfp;
using namespace pre_matfp;

int main(int argc, char *argv[]) {
  int log_level = 1;  // debug
  std::string log_filename;
  int num_subdivide;
  Args args;
  bool is_shown = true;

  CLI::App app{"matfp_pre"};
  app.add_option("input,--input", args.input_surface_path,
                 "Input surface mesh INPUT in .off/.obj/.stl/.ply format. "
                 "(string, required)")
      ->required();
  app.add_option("--sub", args.num_subdivide,
                 "Number of times for subdivision (unsigned int).");
  app.add_option("--log", log_filename, "Log info to given file.");
  app.add_option("--concave", args.thres_concave,
                 "Threshold for detecting concave edges, default 0.18, smaller "
                 "more sensative (double)");
  app.add_option("--convex", args.thres_convex,
                 "Threshold for detecting sharp edges, default 30, smaller "
                 "more sensative (double)");
  app.add_option("--save", args.is_save_model,
                 "Flag for saving preprocessed model (boolean).");
  app.add_option("--scaled", args.is_save_scaled_input,
                 "Flag for saving scaled input model (boolean).");
  app.add_option("--show", is_shown,
                 "Flag for showing polyscope GUI (boolean).");
  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  matfp::Logger::init("pre_matfp", true, log_filename);
  log_level = std::max(0, std::min(6, log_level));
  spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
  spdlog::flush_every(std::chrono::seconds(3));

  // this will make sure all geogram algoriths can be used
  GEO::initialize();
  GEO::CmdLine::import_arg_group("algo");

  PreGuiWindow gui(args);
  if (is_shown) gui.show();

  // will shutdown when program terminates
  // GEO::terminate();
  return 0;
}