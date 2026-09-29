// File: mirrow_boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 15.01.10.
// Description: Program for calculation of rarefaid flows.

#include "mirror_boundary.h"

#include <cmath>
#include <limits>

namespace mc3d {
MirrorBoundary::MirrorBoundary() {}

MirrorBoundary::MirrorBoundary(Point pstn, Point nrml) {
  if ((abs(pstn.x) >= abs(pstn.y)) && (abs(pstn.x) >= abs(pstn.z))) {
    pstn.y = 0;
    pstn.z = 0;
    nrml.y = 0;
    nrml.z = 0;
    if (nrml.x == 0) exit(-1);
  } else if ((abs(pstn.y) >= abs(pstn.x)) && (abs(pstn.y) >= abs(pstn.z))) {
    pstn.x = 0;
    pstn.z = 0;
    nrml.x = 0;
    nrml.z = 0;
    if (nrml.y == 0) exit(-1);
  } else if ((abs(pstn.z) >= abs(pstn.x)) && (abs(pstn.z) >= abs(pstn.y))) {
    pstn.x = 0;
    pstn.y = 0;
    nrml.x = 0;
    nrml.y = 0;
    if (nrml.z == 0) exit(-1);
  }
  nrml.Normalize();
  this->nrml = nrml;
  this->pstn = pstn;
}

MirrorBoundary::~MirrorBoundary() {}

int MirrorBoundary::BoundaryCondition(std::vector<Particle>& particles,
                                      double dt) {
  double t;
  auto data = particles.begin();
  const double Pi = 3.1415926535;
  while (data != particles.end()) {
    t = (pstn - data->position) * nrml;

    if (t < 0 || (t == 0 && data->velocity * nrml > 0)) {
      t /= nrml * nrml;

      data->position += nrml * t * 2.;
      if (!(nrml * Point(0, 1, 1))) {
        data->velocity.x = -data->velocity.x;
      } else if (!(nrml * Point(1, 0, 1))) {
        data->velocity.y = -data->velocity.y;
      } else if (!(nrml * Point(1, 1, 0))) {
        data->velocity.z = -data->velocity.z;
      } else
        return -1;
    }

    // Keep exact upper-face hits inside the half-open cell domain.
    if (t == 0) {
      const double inward = -std::numeric_limits<double>::infinity();
      if (nrml.x > 0) data->position.x = std::nextafter(data->position.x, inward);
      if (nrml.y > 0) data->position.y = std::nextafter(data->position.y, inward);
      if (nrml.z > 0) data->position.z = std::nextafter(data->position.z, inward);
    }
    ++data;
  }
  return 1;
}
}  // namespace mc3d