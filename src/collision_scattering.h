#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>

#include "point.h"
namespace mc3d::detail {
// Uniform solid angle: the cosine, rather than the polar angle, is uniform.
inline Point IsotropicRelativeVelocity(double g, double u, double v) {
  const double mu = 2 * u - 1;
  const double phi = 2 * std::numbers::pi * v;
  const double s = std::sqrt(std::max(0., 1 - mu * mu));
  return {g * s * std::cos(phi), g * s * std::sin(phi), g * mu};
}
}  // namespace mc3d::detail
