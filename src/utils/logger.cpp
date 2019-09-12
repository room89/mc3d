#include "logger.hpp"

namespace logging {

namespace {
double TimeFromStart() {
  using namespace std::chrono;
  static auto start = high_resolution_clock::now();
  auto now = high_resolution_clock::now();
  return duration_cast<milliseconds>(now - start).count() / 1000.0;
}
}  // namespace

Logger::Logger(const std::string& file, const size_t line_number,
               const std::string& level) {
  auto now = std::chrono::system_clock::now();
  std::time_t now_c =
      std::chrono::system_clock::to_time_t(now - std::chrono::hours(24));
  common_stream_ << "tskv\tlevel=" << level << "\tfile=" << file
                 << "(line:" << line_number
                 << ")\ttime=" << std::put_time(std::localtime(&now_c), "%F %T")
                 << "\tfrom_start=" << TimeFromStart();
}

Logger::~Logger() {
  std::cerr << common_stream_.str() << "\ttext=" << text_stream_.str()
            << std::endl;
}

Logger& Logger::operator<<(const LogExtra& extra) {
  for (auto& p : extra) {
    common_stream_ << "\t" << p.first << "=" << p.second;
  }
  return *this;
}

}  // namespace logging
