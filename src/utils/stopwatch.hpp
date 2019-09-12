#include <chrono>
#include <string>

namespace utils {

class Stopwatch {
 public:
  Stopwatch(const std::string& name);
  ~Stopwatch();

  using Clock = std::chrono::high_resolution_clock;
  using TimePoint = Clock::time_point;

 private:
  const std::string name_;

  const TimePoint start_time_;
};

class ScopeWather {
 public:
  ScopeWather(Stopwatch& stopwatch, const std::string& scope_name);
};

}  // namespace utils
