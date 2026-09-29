#include "free_boundary.h"

#include <algorithm>
#include <cmath>
#include <utils/logger.hpp>

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

FreeBoundary::FreeBoundary() {}

FreeBoundary::~FreeBoundary() {}

FreeBoundary::FreeBoundary(Point pstn, Point nrml, unsigned int np, double S,
                           double T, double alpha) {
  if ((abs(pstn.x) >= abs(pstn.y)) && (abs(pstn.x) >= abs(pstn.z))) {
    pstn.y = 0;
    pstn.z = 0;
    nrml.y = 0;
    nrml.z = 0;
    if (nrml.x == 0)
      exit(unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR);
  } else if ((abs(pstn.y) >= abs(pstn.x)) && (abs(pstn.y) >= abs(pstn.z))) {
    pstn.x = 0;
    pstn.z = 0;
    nrml.x = 0;
    nrml.z = 0;
    if (nrml.y == 0)
      exit(unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR);
  } else if ((abs(pstn.z) >= abs(pstn.x)) && (abs(pstn.z) >= abs(pstn.y))) {
    pstn.x = 0;
    pstn.y = 0;
    nrml.x = 0;
    nrml.y = 0;
    if (nrml.z == 0)
      exit(unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR);
  }
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
  const double incoming =
      dt * np * sqrt(T / (Pi * 2)) *
      (exp(-Vn * Vn / (2 * T)) +
       sqrt(Pi) * (Vn / sqrt(2 * T)) * (1 + std::erf(Vn / sqrt(2 * T)))) /
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
