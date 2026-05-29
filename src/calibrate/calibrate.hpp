#pragma once
#include <string>

namespace genopack::calibrate {

int run_calibrate(const std::string& archive,
                  const std::string& checkm2,
                  const std::string& output,
                  int threads);

} // namespace genopack::calibrate
