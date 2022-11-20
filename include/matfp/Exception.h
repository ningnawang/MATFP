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
