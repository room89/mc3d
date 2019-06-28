#pragma once

#include "free_boundary.h"

namespace mc3d {
class giper_free_boundary : public free_boundary {
 public:
  giper_free_boundary(void);
  giper_free_boundary(point pstn, point nrml, unsigned int np, double S,
                      double T, double alpha);
  ~giper_free_boundary(void);
  int bondary_condition(deque<particle>* cluster_particle, double dt);
};
}  // namespace mc3d
