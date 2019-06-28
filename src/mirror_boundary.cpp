// File: mirrow_boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 15.01.10.
// Description: Program for calculation of rarefaid flows.

#include "mirror_boundary.h"

namespace mc3d {
mirror_boundary::mirror_boundary(void) {}

mirror_boundary::mirror_boundary(point pstn, point nrml) {
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
  nrml.nrmlz();
  this->nrml = nrml;
  this->pstn = pstn;
}

mirror_boundary::~mirror_boundary(void) {}

int mirror_boundary::bondary_condition(deque<particle>* particles, double dt) {
  double t;
  deque<particle>::iterator data = particles->begin();
  const double Pi = 3.1415926535;
  while (data != particles->end()) {
    // t = ((pstn.x - data->x) * nrml.x + (pstn.y - data->y) * nrml.y + (pstn.z
    // - data->z) * nrml.z);
    t = (pstn - data->position) * nrml;

    if (t < 0) {
      double rmt = 1 / double(RAND_MAX);
      t /= nrml * nrml;

      /*data->x += 2 * nrml.x * t;
      data->y += 2 * nrml.y * t;
      data->z += 2 * nrml.z * t;*/
      data->position += nrml * t * 2.;
      if (!(nrml * point(0, 1, 1))) {
        data->velocity.x = -data->velocity.x;
      } else if (!(nrml * point(1, 0, 1))) {
        data->velocity.y = -data->velocity.y;
      } else if (!(nrml * point(1, 1, 0))) {
        data->velocity.z = -data->velocity.z;
      } else
        return -1;
      /*double v = sqrt(data->u * data->u + data->v * data->v + data->w *
      data->w); double rn1 = double(std::rand()) * rmt; double rn2 =
      double(std::rand()) * rmt; data->u = v * sin(Pi * rn1) * cos(2 * Pi *
      rn2); data->v = v * sin(Pi * rn1) * sin(2 * Pi * rn2); data->w = v *
      cos(Pi * rn1); point vel(data->u, data->v, data->w); if(vel * nrml > 0)
      {
              data->u = -data->u;
              data->v = -data->v;
              data->w = -data->w;
      }*/
    }

    data++;
  }
  return 1;
}
}  // namespace mc3d