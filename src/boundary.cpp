// File: boundary.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 16.01.09.
// Description: Program for calculation of rarefaid flows.

#include "boundary.h"

namespace mc3d {
void boundary::set_normal(point a) { nrml = a; }

void boundary::set_position(point a) { pstn = a; }
}  // namespace mc3d