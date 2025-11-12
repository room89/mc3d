#include "utils.hpp"

#include <atomic>
#include <random>

namespace {
std::mt19937 CreateEngine(uint32_t seed_base, uint32_t counter) {
  return std::mt19937(seed_base + counter * 9973u);
}
}  // namespace

namespace utils {
namespace {
std::atomic<uint32_t> seed_base{[]() {
  std::random_device rd;
  return rd();
}()};
std::atomic<uint32_t> seed_counter{0};
thread_local std::uniform_real_distribution<double> uniform01_distribution{0.0,
                                                                           1.0};
}  // namespace

std::mt19937& RandomEngine() {
  thread_local std::mt19937 engine =
      CreateEngine(seed_base.load(std::memory_order_relaxed),
                   seed_counter.fetch_add(1, std::memory_order_relaxed));
  return engine;
}

void SeedRandom(uint32_t seed) {
  seed_base.store(seed, std::memory_order_relaxed);
  seed_counter.store(0, std::memory_order_relaxed);
  // Re-seeding existing thread-local engines requires reinitialization.
  // Force recreation on next use by resetting the thread-local counter.
  thread_local bool reset = []() {
    RandomEngine() =
        CreateEngine(seed_base.load(std::memory_order_relaxed),
                     seed_counter.fetch_add(1, std::memory_order_relaxed));
    uniform01_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
    return true;
  }();
  (void)reset;
}

double RandomDouble(double min, double max) {
  thread_local std::uniform_real_distribution<double> distribution;
  if (distribution.a() != min || distribution.b() != max) {
    distribution = std::uniform_real_distribution<double>(min, max);
  }
  return distribution(RandomEngine());
}

double Random01() { return uniform01_distribution(RandomEngine()); }

}  // namespace utils
