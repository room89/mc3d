#include "runtime_control.h"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <thread>
#include <utility>
#include "utils/logger.hpp"

namespace mc3d {
RuntimeControl::RuntimeControl(std::string path, std::size_t capacity)
    : path_(std::move(path)), capacity_(capacity) {
  if (!capacity_) throw std::invalid_argument("Runtime control requires worker capacity");
}

void RuntimeControl::Refresh() {
  if (path_.empty()) return;
  std::ifstream file(path_, std::ios::binary);
  if (!file) return;  // Missing/unreadable file retains the last accepted mode.
  char buffer[65];
  file.read(buffer, sizeof(buffer));
  std::string request(buffer, static_cast<std::size_t>(file.gcount()));
  if (file.bad()) return;
  if (file.gcount() == sizeof(buffer)) request = "<oversized request>";
  const auto first = request.find_first_not_of(" \t\r\n");
  const auto last = request.find_last_not_of(" \t\r\n");
  request = first == std::string::npos ? "" : request.substr(first, last - first + 1);
  if (request == last_request_) return;
  last_request_ = request;
  if (request != "full" && request != "balanced" && request != "quiet" && request != "pause") {
    LOG_WARNING() << "Invalid runtime control mode; retaining " << mode_;
    return;
  }
  if (mode_ == request) return;
  mode_ = request;
  LOG_INFO() << "Runtime load: " << mode_ << "; workers=" << Threads()
             << "; duty=" << DutyCycle();
}

std::size_t RuntimeControl::Threads() const {
  if (mode_ == "quiet") return 1;
  if (mode_ == "balanced") return capacity_ / 2 + capacity_ % 2;
  return capacity_;
}

double RuntimeControl::DutyCycle() const {
  if (mode_ == "quiet") return 0.25;
  if (mode_ == "balanced") return 0.6;
  if (mode_ == "pause") return 0;
  return 1;
}

std::size_t RuntimeControl::WaitBeforeStep(std::chrono::steady_clock::duration previous_work, const std::function<bool()>& should_stop) {
  if (path_.empty()) return capacity_;
  using Clock = std::chrono::steady_clock;
  const auto started = Clock::now();
  for (;;) {
    if (should_stop && should_stop()) return Threads();
    Refresh();
    const auto duty = DutyCycle();
    if (duty > 0) {
      const auto rest = std::chrono::duration<double>(previous_work) * (1 / duty - 1);
      const std::chrono::duration<double> remaining = rest - (Clock::now() - started);
      if (remaining.count() <= 0) return Threads();
      std::this_thread::sleep_for(std::min(remaining, std::chrono::duration<double>(0.2)));
    } else {
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
    }
  }
}
}  // namespace mc3d
