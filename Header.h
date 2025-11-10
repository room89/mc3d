// File: Header.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of freemolecular flows.#pragma once

double F(double x, double y, double z) {
  if (x <= 0.5) return 1.2;
  return 1;
}