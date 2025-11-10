// File: Polygon.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.4.0
// Last modified: 11.10.10.
// Description: Program for calculation of rarefied flows.

// File: Polygon.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of rarefied flows.

#pragma once

#include <utility>

#include "point.h"

namespace mc3d {
struct Polygon {
  Point p1;
  Point p2;
  Point p3;
  Point normal;
  Point force;
  double flux;
  double S;
  // double colision(Particle a, double dt);
  Polygon();
  Polygon(Point p1, Point p2, Point p3, Point nrml);
  Polygon(Point p1, Point p2, Point p3);
  ~Polygon();
  void Set(Point p1, Point p2, Point p3, Point nrml);
  Point GetP1() const;
  Point GetP2() const;
  Point GetP3() const;
  Point GetNormal() const;
  Point GetGmt() const;
  Polygon* GetPtr();
  double GetLmax() const;
  std::pair<Polygon*, Polygon*> Divide();
  void Print();
  void Move(const Point& a);
  void Scale(double e);
  bool Fix();
  double DistanceToPoint(Point p);
};
};  // namespace mc3d
