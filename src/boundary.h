// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <complex>
#include <deque>
#include "exit_code.h"
#include "particle.h"
#include "point.h"

namespace mc3d {
class boundary {
 protected:
  // char type[6];
  point nrml;  //направление граничного условия
  point pstn;  //положение граничного условия
 public:
  virtual int bondary_condition(deque<particle>* a, double dt = 0) = 0;
  virtual ~boundary() {}
  void set_position(point a);  //установка места граничных условий
  void set_normal(point a);  //установка направления граничных условий
};
}  // namespace mc3d
