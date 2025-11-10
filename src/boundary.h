// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <vector>

#include "particle.h"
#include "point.h"

namespace mc3d {
class Boundary {
 protected:
  // char type[6];
  Point nrml;  // направление граничного условия
  Point pstn;  // положение граничного условия
 public:
  virtual int BoundaryCondition(std::vector<Particle>& particles,
                                double dt = 0) = 0;
  virtual ~Boundary() {}
  void SetPosition(Point a);  // установка места граничных условий
  void SetNormal(Point a);    // установка направления граничных условий
};
}  // namespace mc3d
