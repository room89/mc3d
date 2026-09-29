#include "utils.hpp"

#include <atomic>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

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

double IncomingFlux(double drift, double temperature) {
  if (!(temperature > 0) || !std::isfinite(temperature) || !std::isfinite(drift))
    throw std::invalid_argument("Inflow requires finite drift and positive temperature");
  const double sigma = std::sqrt(temperature);
  const double s = drift / sigma;
  return std::max(0.0, sigma / std::sqrt(2 * std::numbers::pi) *
                          std::exp(-0.5 * s * s) +
                          0.5 * drift * std::erfc(-s / std::sqrt(2.0)));
}

double IncomingNormalSpeed(double drift, double temperature) {
  if (!(temperature > 0) || !std::isfinite(temperature) || !std::isfinite(drift))
    throw std::invalid_argument("Inflow requires finite drift and positive temperature");
  const double sigma = std::sqrt(temperature);
  auto log_uniform = [] {
    return std::log(std::max(Random01(), std::numeric_limits<double>::min()));
  };
  if (drift <= 0) {
    for (;;) {
      if (drift < -sigma) {
        // Gamma(shape=2, scale=T/|drift|) envelope, efficient for outflow tails.
        const double w = -(temperature / -drift) * (log_uniform() + log_uniform());
        if (w > 0 && log_uniform() <= -w * w / (2 * temperature)) return w;
      } else {
        // Rayleigh envelope; at zero drift this is the exact distribution.
        const double w = sigma * std::sqrt(-2 * log_uniform());
        if (w > 0 && log_uniform() <= drift * w / temperature) return w;
      }
    }
  }
  // A Gaussian of variance 2*T bounds the flux-weighted Maxwellian.
  // The envelope maximum solves w*(w-drift)=2*T. No velocity cutoff is used.
  const double mode = 0.5 * (drift + std::hypot(drift, std::sqrt(8 * temperature)));
  const double shift = mode - drift;
  std::normal_distribution<double> proposal(drift, std::sqrt(2 * temperature));
  for (;;) {
    const double w = proposal(RandomEngine());
    if (w <= 0) continue;
    const double delta = w - drift;
    const double log_accept = std::log(w / mode) -
                             (delta * delta - shift * shift) / (4 * temperature);
    if (log_uniform() <= log_accept) return w;
  }
}

}  // namespace utils
