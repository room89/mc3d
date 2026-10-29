#include "utils.hpp"

#include <atomic>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <random>
#include <stdexcept>

namespace utils {
double IncomingFlux(double drift, double temperature) {
  if (!(temperature > 0) || !std::isfinite(temperature) || !std::isfinite(drift))
    throw std::invalid_argument("Inflow requires finite drift and positive temperature");
  const double sigma = std::sqrt(temperature);
  const double s = drift / sigma;
  return std::max(0.0, sigma / std::sqrt(2 * std::numbers::pi) *
                          std::exp(-0.5 * s * s) +
                          0.5 * drift * std::erfc(-s / std::sqrt(2.0)));
}

double IncomingNormalSpeed(double drift, double temperature, RandomStream& random) {
  if (!(temperature > 0) || !std::isfinite(temperature) || !std::isfinite(drift))
    throw std::invalid_argument("Inflow requires finite drift and positive temperature");
  const double sigma = std::sqrt(temperature);
  auto log_uniform = [&random] {
    return std::log(std::max(random.Uniform(), std::numeric_limits<double>::min()));
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
    const double w = proposal(random.Engine());
    if (w <= 0) continue;
    const double delta = w - drift;
    const double log_accept = std::log(w / mode) -
                             (delta * delta - shift * shift) / (4 * temperature);
    if (log_uniform() <= log_accept) return w;
  }
}

}  // namespace utils
