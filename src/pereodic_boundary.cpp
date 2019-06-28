// File: pereodic_boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#include "pereodic_boundary.h"

namespace mc3d {
pereodic_boundary::pereodic_boundary(void) {}

pereodic_boundary::~pereodic_boundary(void) {}

int pereodic_boundary::bondary_condition(deque<particle>* particles,
                                         double dt) {
  deque<particle>::iterator data = particles->begin();
  while (data != particles->end()) {
    while ((data->position - pstn) * nrml > 0) {
      data->position += mixing;
    }
    data++;
  }
  return 0;
}

pereodic_boundary::pereodic_boundary(point pstn, point nrml, point mixing) {
  this->nrml = nrml;
  this->pstn = pstn;
  this->mixing = mixing;
}

void pereodic_boundary::set_mixing(point b) { mixing = b; }
}  // namespace mc3d