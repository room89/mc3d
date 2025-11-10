#pragma once

#include <cstdint>
#include <random>

namespace utils {

std::mt19937& RandomEngine();

void SeedRandom(uint32_t seed);

double RandomDouble(double min, double max);

double Random01();

}  // namespace utils
