#include "utils.hpp"

namespace utils {
double RandomDouble(double max, double min) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(min, max);
  return dis(gen);
}
}  // namespace utils
