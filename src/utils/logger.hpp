#pragma once

#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>

namespace logging {

using LogExtra = std::map<std::string, std::string>;

class Logger {
 public:
  Logger(const std::string& file, const size_t line_number,
         const std::string& level);
  ~Logger();
  Logger& operator<<(const LogExtra& extra);

  template <typename T>
  Logger& operator<<(const T& extra) {
    text_stream_ << extra;
    return *this;
  }

 private:
  std::stringstream common_stream_;
  std::stringstream text_stream_;
};

}  // namespace logging

#define LOG_INFO() (logging::Logger(__FILE__, __LINE__, "INFO"))
#define LOG_ERROR() (logging::Logger(__FILE__, __LINE__, "ERROR"))
#define LOG_WARNING() (logging::Logger(__FILE__, __LINE__, "WARNING"))
#define LOG_DEBUG() (logging::Logger(__FILE__, __LINE__, "DEBUG"))
