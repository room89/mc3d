// File: Particle.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include "point.h"

namespace mc3d {
struct Particle {
 public:
  Point velocity;
  Point position;

  Particle();
  Particle(double x, double y, double z, double u, double v, double w);
  Particle(Point position, Point velocity);
  ~Particle();

  void Move(double dt);

  Point GetVelocity() const { return Point(velocity); }
  Point GetPosition() const { return Point(position); }

  double GetU() const;
  double GetV() const;
  double GetW() const;

  double GetX() const;
  double GetY() const;
  double GetZ() const;

  void SetVelocity(Point vel);
  void SetPosition(Point pos);
  void SetX(double x);
  void SetY(double y);
  void SetZ(double z);
  void SetU(double u);
  void SetV(double v);
  void SetW(double w);
};

bool Collision(Particle& a, Particle& b, double& g_max_ch,
               double& frequency_t_ch, double factor);
}  // namespace mc3d