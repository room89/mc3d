// File: particle.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#pragma once

#include "point.h"

namespace mc3d {
struct particle {
 public:
  // double x, y, z;
  // double u, v, w;
  point velocity;
  point position;
  // double r;
  // friend void collision(particle a, particle b);
 public:
  particle(void);
  ~particle(void);
  particle(double x, double y, double z, double u, double v, double w);
  particle(point position, point velocity);
  // void print();
  void move(double dt);
  point get_velocity() const { return point(velocity); };
  point get_position() const { return point(position); };
  double get_u();
  double get_v();
  double get_w();
  double get_x();
  double get_y();
  double get_z();
  void set_velocity(point vel);
  void set_position(point pos);
  void set_x(double x);
  void set_y(double y);
  void set_z(double z);
  void set_u(double u);
  void set_v(double v);
  void set_w(double w);
  // friend void collision(particle a, particle b);
};

bool collision(particle* a, particle* b, double& g_max_ch,
               double& frequency_t_ch, double factor);
}  // namespace mc3d