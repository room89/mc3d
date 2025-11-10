// File: var.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 21.12.09.
// Description: Program for calculation of ra refaid.

#pragma once

#include <fstream>
#include <iostream>

#include "cell.h"
#include "point.h"

namespace mc3d {
struct Var {
  Point apex;
  Point size;
  Point N;
  Point velocity;
  double S;
  double alpha;
  double T;

 public:
  friend std::ostream& operator<<(std::ostream&, const Var&);
  Cell* MakeCell();
};
}  // namespace mc3d
