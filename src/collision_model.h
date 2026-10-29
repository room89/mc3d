#pragma once
#include <vector>
#include "particle.h"
#include "utils/utils.hpp"
namespace mc3d {
struct CollisionModel {
  static void Apply(std::vector<Particle>& particles, Point velocity, double length, double local_kn, double dt, utils::RandomStream& random);
};
}
