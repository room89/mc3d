#include "free_boundary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utils/logger.hpp>
#include <utils/utils.hpp>

namespace mc3d {

FreeBoundary::FreeBoundary() {}

FreeBoundary::~FreeBoundary() {}

FreeBoundary::FreeBoundary(Point pstn, Point nrml, unsigned int np, double S,
                           double T, double alpha) {
  // A plane through the origin still has an axis: use its normal, not position.
  if ((nrml.x != 0) + (nrml.y != 0) + (nrml.z != 0) != 1 ||
      !std::isfinite(nrml.Mod()) || nrml.Mod() == 0)
    throw std::invalid_argument("Open boundary requires an axis-aligned normal");
  nrml.Normalize();
  this->nrml = nrml;
  this->pstn = pstn;

  this->np = np;
  this->T = T;
  this->V.Set(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
  this->S = S;

  Vn = -(V * nrml);

  LOG_DEBUG() << "Velocity: (" << V.x << ", " << V.y << ", " << V.z
              << ") Vn: " << Vn << " pstn " << pstn;
}

void FreeBoundary::AddCell(std::deque<Cell>& cluster_cells) {
  for (auto& cell : cluster_cells) {
    if (abs(nrml * (cell.GetCenter() - pstn)) <
        0.6 * abs(nrml * cell.GetSize())) {
      cells_.push_back({cell, 0});
    }
  }
}

unsigned int FreeBoundary::AccumulateInflow(InflowCell& entry, double dt,
                                          unsigned int minimum_batch) {
  if (dt <= 0 || np == 0) return 0;
  const double incoming = dt * np * utils::IncomingFlux(Vn, T) /
                          abs(entry.cell.get().GetSize() * nrml);
  // Keep the fractional flux in this cell across variable time steps.
  entry.pending_particles += std::max(0.0, incoming);
  if (entry.pending_particles < minimum_batch) return 0;
  const auto count = static_cast<unsigned int>(entry.pending_particles);
  entry.pending_particles -= count;
  return count;
}

int FreeBoundary::BoundaryCondition(std::vector<Particle>& /*cluster_particle*/,
                                    double dt) {
  for (auto& entry : cells_) {
    // GenerateFreeRandom normalizes the batch temperature, requiring at least
    // three particles. Small inflows are delayed, never discarded.
    const auto count = AccumulateInflow(entry, dt, 3);
    if (count > 0) entry.cell.get().GenerateFreeRandom(count, T, V, nrml);
  }

  return 1;
}

void FreeBoundary::SetNp(unsigned int np) { this->np = np; }

}  // namespace mc3d
