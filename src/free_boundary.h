// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 16.12.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include "boundary.h"
#include "cell.h"

namespace mc3d {
class FreeBoundary : public Boundary {
 protected:
  deque<mc3d::Cell*> cells_ptr;
  deque<mc3d::Cell*>::iterator cell_iter;
  unsigned int np;
  double Vn;
  double T;
  double S;
  Point V;

 public:
  virtual void AddCell(deque<Cell>& cluster_cells);
  virtual void SetNp(unsigned int np);
  int BoundaryCondition(deque<Particle>* cluster_particle, double dt);
  FreeBoundary(Point pstn, Point nrml, unsigned int np, double S, double T,
               double alpha = 0);
  FreeBoundary();
  ~FreeBoundary();
};
}  // namespace mc3d
