//
// This file is part of the libWetHair open source project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Copyright 2023 Fredrik Salomonsson

#ifndef LIBWETHAIR_CORE_LOG_UTILITIES_H_
#define LIBWETHAIR_CORE_LOG_UTILITIES_H_

#include <filesystem>
#include <string_view>

namespace libwethair {

/// Generate a path to a file in the temp directory.
///
/// The path will be <temp directory>/<prefix><timestamp>.log
std::filesystem::path tempLogFile(const std::string_view prefix,
                                  const char* timestamp_format = "%Y%m%dT%H%M%S");

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_LOG_UTILITIES_H_

