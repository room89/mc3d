#include "giper_free_boundary.h"

namespace mc3d {

HyperFreeBoundary::HyperFreeBoundary() {}

HyperFreeBoundary::HyperFreeBoundary(Point pstn, Point nrml, unsigned int np,
                                     double S, double T, double alpha)
    : FreeBoundary(pstn, nrml, np, S, T, alpha) {}

HyperFreeBoundary::~HyperFreeBoundary() {}

int HyperFreeBoundary::BoundaryCondition(
    std::vector<Particle>& /*cluster_particle*/, double dt) {
  for (auto& entry : cells_) {
    const auto count = AccumulateInflow(entry, dt, 1);
    if (count > 0) entry.cell.get().GenerateHyperFreeRandom(count, V, T);
  }

  return 1;
}
}  // namespace mc3d
