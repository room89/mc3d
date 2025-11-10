// File: var.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 22.12.09.
// Description: Program for calculation of rarefaid flows.

#include "var.h"

#include <iostream>

namespace mc3d {
Cell* Var::MakeCell() {
  Cell* new_cell = new Cell;

  new_cell->SetTemperature(T);
  new_cell->SetVelocity(velocity);
  return nullptr;
}

std::ostream& operator<<(std::ostream& o, const Var& c) {
  o << c.apex + c.size / 2 << c.N << "\t" << c.T << c.velocity << "\t";
  return o;
}
}  // namespace mc3d