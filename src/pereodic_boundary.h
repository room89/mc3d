// File: pereodic_boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 55.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once
#include "boundary.h"

namespace mc3d {
class pereodic_boundary : public boundary {
 private:
  point mixing;

 public:
  pereodic_boundary(point pstn, point nrml, point mixing);
  void set_mixing(point b);
  int bondary_condition(deque<particle>* cluster_particle, double dt);
  pereodic_boundary(void);
  ~pereodic_boundary(void);
};
}  // namespace mc3d