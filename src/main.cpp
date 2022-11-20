// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOBJ.h>
#include <igl/write_triangle_mesh.h>
#include <matfp/Common.h>
#include <matfp/DisableWarnings.h>
#include <matfp/EnableWarnings.h>
#include <matfp/Logger.h>
#include <matfp/ThreeDimensionalShape.h>

#include <CLI/CLI.hpp>

#include "matfp/GuiWindow.h"

using namespace matfp;
int main(int argc, char *argv[]) {
  int log_level = 1;  // debug
  std::string log_filename;
  std::string mat_path;
  Args args;

  CLI::App app{"matfp"};
  // More details please see Args.h
  app.add_option("input,--input", args.input_surface_path,
                 "Input surface mesh INPUT in .geogram format. "
                 "(string, required)")
      ->required();
  app.add_option("mat,--ma", mat_path, "MAT surface in .ma format. (string)");
  app.add_option(
      "downsample,--ds", args.downsample_percentage,
      "Downsample percentage when input triangles are dense. Once set, "
      "we will not use rsample. Larger the denser. (double)");
  app.add_option("rsample,--r", args.rsample,
                 "Control random sampling rate, r-sample from Nina Amenta "
                 "(check power crust). Smaller the denser. (double)");
  app.add_option(
      "cc_len_eps,--cc", args.cc_len_eps,
      "To specify pin points on concave edge with length cc_len_eps. "
      "Default=0.03. Performs better when cc_len_eps is smaller but needs more "
      "processing time because more spheres inserted to concave lines. Scaled "
      "in [0,10]. See Args.h for more detail. "
      "(double)");
  app.add_option("--level", log_level,
                 "Log level (0 = most verbose, 6 = off).");

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    return app.exit(e);
  }

  matfp::Logger::init("matfp", true, log_filename);
  log_level = std::max(0, std::min(6, log_level));
  spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
  spdlog::flush_every(std::chrono::seconds(3));

  // this will make sure all geogram algoriths can be used
  GEO::initialize();
  GEO::CmdLine::import_arg_group("algo");

  // 3D
  GuiWindow gui(args);
  if (!mat_path.empty()) {
    gui.set_mat_path(mat_path);
  }
  // gui.automation();
  gui.show();

  // will shutdown when program terminates
  // GEO::terminate();
  return 0;
}