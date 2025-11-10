// File: var.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.3.3
// Last modified: 9.07.09.
// Description: Program for calculation of freemolecular flows.

#pragma once

#include <iostream>
using namespace std;

struct var {
  double u;
  double v;
  double w;
  double t;
  unsigned int n;

 public:
  void set(double u, double v, double w, double t, unsigned int n);
  // operator<<;
  var(void);
  var(double u, double v, double w, double t, unsigned int n);
  ~var(void);
  friend ostream &operator<<(ostream &o, const var &c);
};
