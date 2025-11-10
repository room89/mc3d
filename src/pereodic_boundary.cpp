// File: pereodic_boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#include "pereodic_boundary.h"

namespace mc3d {
PeriodicBoundary::PeriodicBoundary() {}

PeriodicBoundary::~PeriodicBoundary() {}

int PeriodicBoundary::BoundaryCondition(deque<Particle>* particles, double dt) {
  deque<Particle>::iterator data = particles->begin();
  while (data != particles->end()) {
    while ((data->position - pstn) * nrml > 0) {
      data->position += mixing;
    }
    data++;
  }
  return 0;
}

PeriodicBoundary::PeriodicBoundary(Point pstn, Point nrml, Point mixing) {
  this->nrml = nrml;
  this->pstn = pstn;
  this->mixing = mixing;
}

void PeriodicBoundary::SetMixing(Point b) { mixing = b; }
}  // namespace mc3d