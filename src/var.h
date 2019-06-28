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

using namespace std;

namespace mc3d {
struct var {
  point apex;
  point size;
  point N;
  point velocity;
  double S;
  double alpha;
  double T;

 public:
  friend ostream& operator<<(ostream&, var);
  cell* make_cell();
};
}  // namespace mc3d
