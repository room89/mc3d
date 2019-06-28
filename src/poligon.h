// File: poligon.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.4.0
// Last modified: 11.10.10.
// Description: Program for calculation of rarefied flows.

// File: poligon.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of rarefied flows.

#pragma once

#include <utility>
#include "exception.h"
#include "particle.h"
#include "point.h"

namespace mc3d {
struct poligon {
  point p1;
  point p2;
  point p3;
  point normal;
  point force;
  double flux;
  double S;
  // double colision(particle a, double dt);
  poligon(void);
  poligon(point p1, point p2, point p3, point nrml);
  poligon(point p1, point p2, point p3);
  ~poligon(void);
  void set(point p1, point p2, point p3, point nrml);
  point get_p1() const;
  point get_p2() const;
  point get_p3();
  point get_normal();
  point get_gmt() const;
  poligon* get_ptr();
  double get_Lmax();
  std::pair<poligon*, poligon*> divide();
  void print();
  void move(point a);
  void scaling(double e);
  bool fix();
  double dist_to_point(point p);
};
};  // namespace mc3d
