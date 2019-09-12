#include "giper_free_boundary.h"

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

giper_free_boundary::giper_free_boundary(void) {}

giper_free_boundary::giper_free_boundary(point pstn, point nrml,
                                         unsigned int np, double S, double T,
                                         double alpha)
    : FreeBoundary(pstn, nrml, np, S, T, alpha) {}

giper_free_boundary::~giper_free_boundary(void) {}

int giper_free_boundary::bondary_condition(
    std::deque<particle>* cluster_particle, double dt) {
  if (Vn < 0) return 0;

  std::deque<cell*>::iterator cell_iter = this->cells_ptr.begin();

  unsigned int N;
  while (cell_iter != cells_ptr.end()) {
    point cell_size = (*cell_iter)->get_size();

    N = static_cast<unsigned int>(
        dt * np * sqrt(T / (Pi * 2)) *
        (exp(-Vn * Vn / (2 * T)) +
         sqrt(Pi) * (Vn / sqrt(2 * T)) * (1 + erf(Vn / sqrt(2 * T)))) /
        abs(cell_size * nrml));

    if (N > 2) (*cell_iter)->generate_giper_free_random(N, V, T);

    cell_iter++;
  }

  return 1;
}
}  // namespace mc3d
