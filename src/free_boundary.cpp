#include "free_boundary.h"

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
      cells_.emplace_back(cell);
    }
  }
}

int FreeBoundary::BoundaryCondition(std::vector<Particle>& /*cluster_particle*/,
                                    double dt) {
  unsigned int N;
  for (Cell& cell : cells_) {
    Point cell_size = cell.GetSize();

    N = static_cast<unsigned int>(
        dt * np * sqrt(T / (Pi * 2)) *
        (exp(-Vn * Vn / (2 * T)) +
         sqrt(Pi) * (Vn / sqrt(2 * T)) * (1 + std::erf(Vn / sqrt(2 * T)))) /
        abs(cell_size * nrml));

    if (N > 2) cell.GenerateFreeRandom(N, T, V, nrml);
  }

  return 1;
}

void FreeBoundary::SetNp(unsigned int np) { this->np = np; }

}  // namespace mc3d
