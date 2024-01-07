//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2023 Fredrik Salomonsson

#include "LogUtilities.h"

#include <chrono>
#include <filesystem>
#include <sstream>
#include <string_view>

namespace libwethair {

std::filesystem::path tempLogFile(const std::string_view prefix,
                                  const char* timestamp_format) {
  const auto now =
      std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  std::stringstream name;
  name << prefix << std::put_time(std::localtime(&now), timestamp_format)
       << ".log";
  return std::filesystem::temp_directory_path() / name.str();
}
}  // namespace libwethair
