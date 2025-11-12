// File: pereodic_boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#include "pereodic_boundary.h"

namespace mc3d {
PeriodicBoundary::PeriodicBoundary() {}

PeriodicBoundary::~PeriodicBoundary() {}

int PeriodicBoundary::BoundaryCondition(std::vector<Particle>& particles,
                                        double dt) {
  auto data = particles.begin();
  while (data != particles.end()) {
    while ((data->position - pstn) * nrml > 0) {
      data->position += mixing;
    }
    ++data;
  }
  return 0;
}

PeriodicBoundary::PeriodicBoundary(Point pstn, Point nrml, Point mixing) {
  this->nrml = nrml;
  this->pstn = pstn;
  this->mixing = mixing;
}

void PeriodicBoundary::SetMixing(Point b) { mixing = b; }

const Point& PeriodicBoundary::GetMixing() const { return mixing; }
}  // namespace mc3d