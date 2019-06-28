// File: var.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 22.12.09.
// Description: Program for calculation of rarefaid flows.

#include "var.h"
#include <fstream>
#include <iostream>

namespace mc3d {
cell* var::make_cell() {
  cell* new_cell = new cell;

  new_cell->set_t(T);
  new_cell->set_vel(velocity);
  return NULL;
}

std::ostream& operator<<(std::ostream& o, const var c) {
  o << c.apex + c.size / 2 << c.N << "\t" << c.T << c.velocity << "\t";
  return o;
}
}  // namespace mc3d