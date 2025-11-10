// File: mirror_boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan

// Last modified: 26.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once
#include "boundary.h"
#include "point.h"

namespace mc3d {
class MirrorBoundary : public Boundary {
 public:
  MirrorBoundary();
  MirrorBoundary(Point pstn, Point nrml);
  ~MirrorBoundary();
  int BoundaryCondition(deque<Particle>* cluster_particle, double dt);
};
};  // namespace mc3d