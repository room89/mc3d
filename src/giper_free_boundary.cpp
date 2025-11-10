#include "giper_free_boundary.h"

#include <cmath>

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

HyperFreeBoundary::HyperFreeBoundary() {}

HyperFreeBoundary::HyperFreeBoundary(Point pstn, Point nrml, unsigned int np,
                                     double S, double T, double alpha)
    : FreeBoundary(pstn, nrml, np, S, T, alpha) {}

HyperFreeBoundary::~HyperFreeBoundary() {}

int HyperFreeBoundary::BoundaryCondition(
    std::deque<Particle>& /*cluster_particle*/, double dt) {
  if (Vn < 0) return 0;

  unsigned int N;
  for (Cell& cell : cells_) {
    Point cell_size = cell.GetSize();

    N = static_cast<unsigned int>(
        dt * np * sqrt(T / (Pi * 2)) *
        (exp(-Vn * Vn / (2 * T)) +
         sqrt(Pi) * (Vn / sqrt(2 * T)) * (1 + std::erf(Vn / sqrt(2 * T)))) /
        abs(cell_size * nrml));

    if (N > 2) cell.GenerateHyperFreeRandom(N, V, T);
  }

  return 1;
}
}  // namespace mc3d
