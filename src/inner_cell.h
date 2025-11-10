// File: inner_cell.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.5.1
// Last modified: 18.05.11.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include "cell.h"

namespace mc3d {
class InnerCell : public Cell {
 public:
  InnerCell();
  ~InnerCell();
};
}  // namespace mc3d
