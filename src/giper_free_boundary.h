#pragma once

#include "free_boundary.h"

namespace mc3d {
class HyperFreeBoundary : public FreeBoundary {
 public:
  HyperFreeBoundary();
  HyperFreeBoundary(Point pstn, Point nrml, unsigned int np, double S, double T,
                    double alpha = 0);
  ~HyperFreeBoundary();
  int BoundaryCondition(deque<Particle>& cluster_particle, double dt);
};
}  // namespace mc3d
