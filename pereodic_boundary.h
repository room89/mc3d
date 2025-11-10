#pragma once
#include "boundary.h"

class pereodic_boundary : public boundary {
 private:
  point mixing;

 public:
  pereodic_boundary(point pstn, point nrml, point mixing);
  void set_mixing(point b);
  int bondary_condition(particle_array *cell_particle);
  pereodic_boundary(void);
  ~pereodic_boundary(void);
};
