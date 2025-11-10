// File: var.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.1
// Last modified: 22.12.09.
// Description: Program for calculation of rarefaid flows.

#include "var.h"

#include <iostream>

namespace mc3d {
std::unique_ptr<Cell> Var::MakeCell() {
  auto cell = std::make_unique<Cell>();
  cell->SetTemperature(T);
  cell->SetVelocity(velocity);
  return cell;
}

std::ostream& operator<<(std::ostream& o, const Var& c) {
  o << c.apex + c.size / 2 << c.N << "\t" << c.T << c.velocity << "\t";
  return o;
}
}  // namespace mc3d