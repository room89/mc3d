// File: boundary.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 25.10.10.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include <deque>
#include <functional>
#include <limits>
#include <vector>

#include "particle.h"
#include "point.h"
#include "statistical_moments.h"

namespace mc3d {
class Cell;
struct LocalInflow {
  std::reference_wrapper<Cell> cell;
  std::vector<Particle> particles;
};
struct BoundaryInflow {
  std::vector<LocalInflow> local;
  std::vector<Particle> transit;
};
class Boundary {
 protected:
  // char type[6];
  Point nrml;  // направление граничного условия
  Point pstn;  // положение граничного условия
 public:
  Point step_impulse{0, 0, 0};
  double step_heat = 0;
  virtual void BindCells(std::deque<Cell>& cells) {}
  virtual int ProcessOutgoing(std::vector<Particle>& particles, double dt = 0) { return 0; }
  virtual BoundaryInflow GenerateInflow(double dt) { return {}; }
  virtual double TransportTimeLimit() const { return std::numeric_limits<double>::infinity(); }
  virtual ~Boundary() {}
  void SetPosition(Point a);  // установка места граничных условий
  void SetNormal(Point a);    // установка направления граничных условий
  const Point& GetPosition() const;
  const Point& GetNormal() const;
};
}  // namespace mc3d
