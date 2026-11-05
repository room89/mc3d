#pragma once
#include <cstdint>
#include <vector>

#include "particle.h"
#include "utils/utils.hpp"
namespace mc3d {
struct CollisionStats {
  std::uint64_t accepted_collisions = 0;
};
struct CollisionModel {
  static CollisionStats ApplyWithStats(std::vector<Particle>& particles,
                                       Point velocity, double length,
                                       double local_kn, double dt,
                                       utils::RandomStream& random);
  static void Apply(std::vector<Particle>& particles, Point velocity,
                    double length, double local_kn, double dt,
                    utils::RandomStream& random);
};
}  // namespace mc3d
