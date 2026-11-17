#include "giper_free_boundary.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace mc3d {

HyperFreeBoundary::HyperFreeBoundary() {}

HyperFreeBoundary::HyperFreeBoundary(Point pstn, Point nrml, unsigned int np,
                                     double S, double T, double alpha)
    : FreeBoundary(pstn, nrml, np, S, T, alpha) {}

HyperFreeBoundary::~HyperFreeBoundary() {}

BoundaryInflow HyperFreeBoundary::GenerateInflow(double dt) {
  BoundaryInflow result;
  for (auto& entry : cells_) {
    unsigned int count = 0;
    if (stochastic_rounding_) {
      if (dt <= 0 || np == 0) continue;
      const double expected = dt * np * density_ratio_ * utils::IncomingFlux(Vn, T) /
                              abs(entry.cell.get().GetSize() * nrml);
      if (!std::isfinite(expected) || expected < 0 ||
          expected > std::numeric_limits<unsigned int>::max() - 1.0)
        throw std::overflow_error("Reservoir inflow exceeds batch capacity");
      count = static_cast<unsigned int>(expected + entry.cell.get().SampleInflowRounding());
    } else {
      count = AccumulateInflow(entry, dt, 1);
    }
    if (count > 0)
      entry.cell.get().GenerateHyperFreeRandom(count, V, T, nrml, dt, result.transit);
  }

  return result;
}
}  // namespace mc3d
