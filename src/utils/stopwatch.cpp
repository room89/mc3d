#include "stopwatch.hpp"

#include <utils/logger.hpp>

namespace utils {

Stopwatch::Stopwatch(const std::string& name)
    : name_(name), start_time_(Clock::now()) {}

Stopwatch::~Stopwatch() {
  using namespace std::chrono;
  auto now = Clock::now();
  auto t = duration<double>(now - start_time_).count();
  LOG_INFO() << "stopwatch:" << name_ << ";total_time:" << t;
}

}  // namespace utils
