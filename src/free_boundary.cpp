#include "free_boundary.h"

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

free_boundary::free_boundary(void) {}

free_boundary::~free_boundary(void) {}

free_boundary::free_boundary(point pstn, point nrml, unsigned int np, double S,
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
  nrml.nrmlz();
  this->nrml = nrml;
  this->pstn = pstn;

  this->np = np;
  this->T = T;
  this->V.set(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
  this->S = S;

  Vn = -(V * nrml);
}

void free_boundary::add_cell(std::deque<cell*>* cluster_cells) {
  cell_iter = cluster_cells->begin();

  while (cell_iter != cluster_cells->end()) {
    if (abs(nrml * ((*cell_iter)->get_center() - pstn)) <
        0.6 * abs(nrml * (*cell_iter)->get_size())) {
      this->cells_ptr.push_back((*cell_iter)->get_ptr());
    }

    cell_iter++;
  }
}

int free_boundary::bondary_condition(std::deque<particle>* cluster_particle,
                                     double dt) {
  std::deque<cell*>::iterator cell_iter = this->cells_ptr.begin();

  unsigned int N;
  while (cell_iter != cells_ptr.end()) {
    point cell_size = (*cell_iter)->get_size();

    N = static_cast<unsigned int>(
        dt * np * sqrt(T / (Pi * 2)) *
        (exp(-Vn * Vn / (2 * T)) +
         sqrt(Pi) * (Vn / sqrt(2 * T)) * (1 + erf(Vn / sqrt(2 * T)))) /
        abs(cell_size * nrml));

    if (N > 10000) {
      std::cout << "N inpu: " << N << "   " << std::endl;
    }

    if (N > 2) (*cell_iter)->generate_free_random(N, T, this->V, this->nrml);

    cell_iter++;
  }

  return 1;
}

void free_boundary::set_np(unsigned int np) { this->np = np; }

// double free_boundary::erf(double x0)
//{
//	if(x0 < 0.00001) return 0;
//	double x = 0, I = 0, dx = x0 / 100;
//	while(x <= x0)
//	{
//		I += exp(-(x + dx / 2) * (x + dx / 2)) * dx;
//		x += dx;
//	}
//	return I;
//}

inline double erf(double x) {
  const double Pi = 3.1415926536, Pit = 1. / Pi, Pi2 = Pi * Pi, Pi2t = 1. / Pi2,
               SPi = sqrt(Pi), SPit = 1. / SPi;
  double xx = fabs(x);
  double x2 = x * x;
  double w;

  if (xx <= 0.4) {
    w = 2. * SPit * x *
        (1. + x2 * (-0.33333333 +
                    x2 * (0.1 + x2 * (-0.02380952 + 0.00462963 * x2))));
  } else {
    double t = 1. / (1. + 0.3275911 * xx);
    w = 1. -
        t *
            (0.254829592 +
             t * (-0.284496736 +
                  t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) *
            exp(-x * x);
    x >= 0. ? w = w : w = -w;
  }

  return w;
}
}  // namespace mc3d
