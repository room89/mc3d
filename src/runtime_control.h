#pragma once

#include <chrono>
#include <functional>
#include <cstddef>
#include <string>

namespace mc3d {
// Runtime scheduling only: never changes physical time steps or simulation state.
class RuntimeControl {
 public:
  RuntimeControl(std::string path, std::size_t capacity);
  void Refresh();
  std::size_t WaitBeforeStep(std::chrono::steady_clock::duration previous_work,
      const std::function<bool()>& should_stop = {});
  std::size_t Threads() const;
  double DutyCycle() const;
  const std::string& Mode() const { return mode_; }
 private:
  std::string path_;
  std::size_t capacity_;
  std::string mode_ = "full";
  std::string last_request_;
};
}  // namespace mc3d
