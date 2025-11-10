// File: pereodic_boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 55.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once
#include "boundary.h"

namespace mc3d {
class PeriodicBoundary : public Boundary {
 private:
  Point mixing;

 public:
  PeriodicBoundary(Point pstn, Point nrml, Point mixing);
  void SetMixing(Point b);
  int BoundaryCondition(std::vector<Particle>& cluster_particle, double dt);
  PeriodicBoundary();
  ~PeriodicBoundary();
};
}  // namespace mc3d