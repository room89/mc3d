#include "utils.hpp"

#include <atomic>

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
    return true;
  }();
  (void)reset;
}

double RandomDouble(double min, double max) {
  std::uniform_real_distribution<double> dist(min, max);
  return dist(RandomEngine());
}

double Random01() { return RandomDouble(0.0, 1.0); }

}  // namespace utils
