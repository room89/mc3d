#pragma once

#include <cstdint>
#include <random>

namespace utils {

std::mt19937& RandomEngine();

void SeedRandom(uint32_t seed);

double RandomDouble(double min, double max);

double Random01();

// Moments of a drifting Maxwellian crossing a plane, per unit number density.
double IncomingFlux(double inward_drift, double temperature);
// Density proportional to w * exp(-(w-inward_drift)^2/(2*T)), w > 0.
double IncomingNormalSpeed(double inward_drift, double temperature);

}  // namespace utils
