// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 16.12.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <deque>
#include <functional>
#include <vector>

#include "boundary.h"
#include "cell.h"

namespace mc3d {
class FreeBoundary : public Boundary {
 protected:
  struct InflowCell {
    std::reference_wrapper<Cell> cell;
    double pending_particles = 0;
  };
  std::vector<InflowCell> cells_;
  unsigned int AccumulateInflow(InflowCell& entry, double dt,
                               unsigned int minimum_batch);
  unsigned int np;
  double Vn;
  double T;
  double S;
  Point V;

 public:
  virtual void AddCell(std::deque<Cell>& cluster_cells);
  virtual void SetNp(unsigned int np);
  int BoundaryCondition(std::vector<Particle>& cluster_particle, double dt);
  FreeBoundary(Point pstn, Point nrml, unsigned int np, double S, double T,
               double alpha = 0);
  FreeBoundary();
  ~FreeBoundary();
};
}  // namespace mc3d
