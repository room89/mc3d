#include "collision_scattering.h"
// File: Particle.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#include <cmath>
#include <utils/utils.hpp>

#include "particle.h"

namespace mc3d {

bool Collision(Particle& a, Particle& b, double& g_max, double& frequency_t,
               double factor, utils::RandomStream& random) {
  Point vel_1 = a.GetVelocity();
  Point vel_2 = b.GetVelocity();
  double g = (vel_2 - vel_1).Mod();

  bool rtrn = true;

  if (g > g_max) {
    g_max = g;
    frequency_t = factor / g_max;
    rtrn = false;
  }

  double rr = random.Uniform();

  if (g / g_max > rr) {
    double r1 = random.Uniform();
    double r2 = random.Uniform();

    const Point g1 = detail::IsotropicRelativeVelocity(g, r1, r2);

    a.SetVelocity(0.5 * (vel_1 + vel_2 - g1));

    b.SetVelocity(0.5 * (vel_1 + vel_2 + g1));
  }

  return rtrn;
}

}  // namespace mc3d