// File: mirror_boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan

// Last modified: 26.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once
#include "boundary.h"
#include "point.h"

namespace mc3d {
class mirror_boundary : public boundary {
 public:
  mirror_boundary(void);
  mirror_boundary(point pstn, point nrml);
  ~mirror_boundary(void);
  int bondary_condition(deque<particle>* cluster_particle, double dt);
};
};  // namespace mc3d