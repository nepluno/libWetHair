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
#include <fstream>
#include <iostream>
#include <ostream>
#include <string_view>

namespace libwethair {

/// Generate a path to a file in the temp directory.
///
/// The path will be <temp directory>/<prefix><timestamp>.log
std::filesystem::path tempLogFile(const std::string_view prefix,
                                  const char* timestamp_format = "%Y%m%dT%H%M%S");

/// Call write with either an std::ostream of path or stream.
///
/// Where the interface for write is:
///
///   void write(std::ostream&)
///
/// It will try to open the file pointed to by path and pass that to
/// write. If it succeeds it will write a message to stream containing
/// the path name. If it fails to open path it will pass stream to
/// write instead.
template<typename Writer>
void delegateLogOutput(const std::filesystem::path& path,
                       const Writer& write,
                       std::ostream& stream = std::cerr)
{
  if (std::ofstream log(path); log.is_open()) {
    write(log);
    stream << "Saved log to: " << path <<"\n";
  } else {
    write(stream);
  }
}

}  // namespace libwethair

#endif  // LIBWETHAIR_CORE_LOG_UTILITIES_H_

