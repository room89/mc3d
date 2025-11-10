#include "pereodic_boundary.h"

pereodic_boundary::pereodic_boundary(void) {}

pereodic_boundary::~pereodic_boundary(void) {}

int pereodic_boundary::bondary_condition(particle_array* claster_particle) {
  if (!claster_particle->go_to_start()) return 1;
  p_elem* a;
  int i = 0;
  unsigned int k = claster_particle->calc_N();
  while (true) {
    i++;
    a = claster_particle->get_data();
    while ((a->p.x - pstn.x) * nrml.x + (a->p.y - pstn.y) * nrml.y +
               (a->p.z - pstn.z) * nrml.z >
           0) {
      a->p.x += mixing.x;
      a->p.y += mixing.y;
      a->p.z += mixing.z;
    }
    // claster_particle->calc_N();
    if (!claster_particle->go_forward()) break;
  }
  return 0;
}

pereodic_boundary::pereodic_boundary(point pstn, point nrml, point mixing) {
  this->nrml = nrml;
  this->pstn = pstn;
  this->mixing = mixing;
}

void pereodic_boundary::set_mixing(point b) { mixing = b; }