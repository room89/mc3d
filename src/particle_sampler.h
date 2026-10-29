#pragma once
#include "inner_boundary.h"
#include "utils/utils.hpp"
namespace mc3d {
struct ParticleSampler {
  static double Initial(std::vector<Particle>& particles, Point apex, Point size, size_t N, double T, Point velocity, utils::RandomStream& random);
  static double Reservoir(std::vector<Particle>& particles, Point apex, Point size, unsigned int N, double T, Point V, utils::RandomStream& random);
  static double Free(std::vector<Particle>& particles, Point apex, Point size, unsigned int N, double T, Point V, Point normal, utils::RandomStream& random);
  static void HyperFree(Point apex, Point size, unsigned int N, Point V, double T, Point normal, double dt, InnerBoundary& body, std::vector<Particle>& arrivals, utils::RandomStream& random);
};
}
