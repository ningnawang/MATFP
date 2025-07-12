// This file is part of MATFP, a software for computing medial axis transform
// with feature preservation.
//
// Copyright (C) 2022 Ningna Wang <ningna.wang@utdallas.edu>
//
// This Source Code Form is subject to the terms of the MIT license.
//
//
// Only used in Logger.h
#pragma once

#include <spdlog/fmt/bundled/format.h>
#include <spdlog/fmt/bundled/ranges.h>

#include <array>
#include <vector>

// // Formatter for std::array<T, N>
// template <typename T, std::size_t N>
// struct fmt::formatter<std::array<T, N>> {
//   constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
//   template <typename FormatContext>
//   auto format(const std::array<T, N>& arr, FormatContext& ctx) const {
//     auto out = ctx.out();
//     *out++ = '[';
//     for (std::size_t i = 0; i < N; ++i) {
//       if (i > 0) *out++ = ',', *out++ = ' ';
//       out = fmt::format_to(out, "{}", arr[i]);
//     }
//     *out++ = ']';
//     return out;
//   }
// };

// Formatter for std::vector<bool>
template <>
struct fmt::formatter<std::vector<bool>> {
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }
  template <typename FormatContext>
  auto format(const std::vector<bool>& vec, FormatContext& ctx) const {
    auto out = ctx.out();
    *out++ = '[';
    for (std::size_t i = 0; i < vec.size(); ++i) {
      if (i != 0) *out++ = ',', *out++ = ' ';
      out = fmt::format_to(out, "{}", static_cast<bool>(vec[i]));
    }
    *out++ = ']';
    return out;
  }
};
