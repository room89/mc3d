// File: Header.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <complex>
#include <iostream>

namespace mc3d {
double F(double x, double y, double z) {
  if (x < 0) return 0;
  return 0.1;
}

double F2(double x, double y, double z) {
  if (y * y + x * x < z) return 1.;
  // if(y * y + x * x < 2 * z) return .5;
  return 0;
}

double V(double x, double y, double R) {
  return 0.3 * (x * x + y * y) *
         std::exp(0.5 * (1 - (x * x + y * y) / (R * R))) / R;
}
}  // namespace mc3d