#pragma once
#include <cmath>
#include <initializer_list>
#include <limits>
#include <stdexcept>
namespace mc3d {
enum class CollisionKind { HS, VHS, None };
struct CollisionParameters {
  CollisionKind kind = CollisionKind::HS;
  double mass = NAN, diameter = NAN, reference_temperature = NAN;
  double omega = NAN, temperature_unit = NAN;
  void Validate() const {
    if (kind != CollisionKind::VHS) return;
    for (double x : {mass, diameter, reference_temperature, temperature_unit})
      if (!(x > 0) || !std::isfinite(x))
        throw std::invalid_argument(
            "VHS requires finite positive mass, diameter, reference "
            "temperature and temperature unit");
    if (!std::isfinite(omega) || omega < .5 || omega > 1)
      throw std::invalid_argument("VHS omega must be in [0.5, 1]");
  }
  // SI cross section times relative speed; reduced mass is m/2.
  double PhysicalRate(double g) const {
    constexpr double kb = 1.380649e-23;
    return std::acos(-1.) * diameter * diameter / std::tgamma(2.5 - omega) *
           std::pow(4 * kb * reference_temperature / mass, omega - .5) *
           std::pow(g, 2 - 2 * omega);
  }
  // Effective relative speed normalized by lambda_eq(Tref)/L.
  double Rate(double g) const {
    if (kind == CollisionKind::None) return 0;
    if (kind == CollisionKind::HS || omega == .5) return g;
    const double theta = reference_temperature / temperature_unit;
    const double p = 2 - 2 * omega;
    const double moment = std::pow(4 * theta, p / 2) *
                          std::tgamma((p + 3) / 2) / std::tgamma(1.5);
    return std::sqrt(2.) * std::sqrt(8 * theta / std::acos(-1.)) / moment *
           std::pow(g, p);
  }
};
}  // namespace mc3d
