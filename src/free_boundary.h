// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 16.12.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include "boundary.h"
#include "cell.h"

namespace mc3d {
class free_boundary : public boundary {
 protected:
  deque<mc3d::cell*> cells_ptr;
  deque<mc3d::cell*>::iterator cell_iter;
  unsigned int np;
  double Vn;
  double T;
  double S;
  // static double erf(double x);
  point V;

 public:
  virtual void add_cell(deque<cell*>* cluster_cells);
  virtual void set_np(unsigned int np);
  int bondary_condition(deque<particle>* cluster_particle, double dt);
  free_boundary(point pstn, point nrml, unsigned int np, double S, double T,
                double alpha);
  free_boundary(void);
  ~free_boundary(void);
};

inline double erf(double x);
}  // namespace mc3d