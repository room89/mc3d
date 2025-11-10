#include "stopwatch.hpp"

#include <utils/logger.hpp>

namespace utils {

Stopwatch::Stopwatch(const std::string& name)
    : name_(name), start_time_(Clock::now()) {}

Stopwatch::~Stopwatch() {
  using namespace std::chrono;
  auto now = high_resolution_clock::now();
  auto t = duration_cast<milliseconds>(now - start_time_).count() / 1000.0;
  LOG_INFO() << "stopwatch:" << name_ << ";total_time:" << t;
}

}  // namespace utils
