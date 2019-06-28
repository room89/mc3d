// File: particle.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#include "particle.h"
#include <stdlib.h>
#include <complex>
#include <iostream>

namespace mc3d {
particle::particle(void) {
  // this->r = 0.001;
  // this->u = this->v = this->w = this->x = this->y = this->z = 0;
  position.set(0, 0, 0);
  velocity.set(0, 0, 0);
}

particle::particle(double x, double y, double z, double u, double v, double w) {
  // this->r = 0.001;
  // this->u = u;
  // this->v = v;
  // this->w = w;
  // this->x = x;
  // this->y = y;
  // this->z = z;
  position.set(x, y, z);
  velocity.set(u, v, w);
}

particle::particle(point position, point velocity) {
  this->position = position;
  this->velocity = velocity;
}

particle::~particle(void) {}

/*void particle::print()
{
        std::cout<<" x = "<<x<<" y = "<<y<<" z = "<<z<<" u = "<<u<<" v = "<<v<<"
w = "<<w<<std::endl;
}*/

void particle::move(double dt) {
  position += velocity * dt;

  /*x += dt * u;
  y += dt * v;
  z += dt * w;*/
}
bool collision(particle* a, particle* b, double& g_max, double& frequency_t,
               double factor) {
  double rmt = 1 / RAND_MAX;
  double Pi = 3.1415923565;

  /*double u1 = a->u, v1 = a->v, w1 = a->w;
  double u2 = b->u, v2 = b->v, w2 = b->w;
  double gx = u2 - u1, gy = v2 - v1, gz = w2 - w1;
  double g = sqrt(gx * gx + gy * gy + gz * gz);*/

  point vel_1 = a->get_velocity();
  point vel_2 = b->get_velocity();
  double g = (vel_2 - vel_1).mod();

  bool rtrn = true;

  if (g > g_max) {
    g_max = g;
    frequency_t = factor / g_max;
    rtrn = false;
  }

  double rr = std::rand() * rmt;

  if (g / g_max > rr) {
    double r1 = std::rand() * rmt;
    double r2 = std::rand() * rmt;

    double g1z = g * cos(Pi * r1), g1y = g * sin(Pi * r1) * sin(2. * Pi * r2),
           g1x = g * sin(Pi * r1) * cos(2. * Pi * r2);

    point g1(g * sin(Pi * r1) * cos(2. * Pi * r2),
             g * sin(Pi * r1) * sin(2. * Pi * r2), g * cos(Pi * r1));

    /*a->u = 0.5 * (u1 + u2) - 0.5 * g1x;
    a->v = 0.5 * (v1 + v2) - 0.5 * g1y;
    a->w = 0.5 * (w1 + w2) - 0.5 * g1z;*/

    a->set_velocity(0.5 * (vel_1 + vel_2 - g1));

    /*b->u = 0.5 * (u1 + u2) + 0.5 * g1x;
    b->v = 0.5 * (v1 + v2) + 0.5 * g1y;
    b->w = 0.5 * (w1 + w2) + 0.5 * g1z;*/

    b->set_velocity(0.5 * (vel_1 + vel_2 + g1));
  }

  return rtrn;
}

double particle::get_u() { return velocity.x; }

double particle::get_v() { return velocity.y; }

double particle::get_w() { return velocity.z; }

double particle::get_x() { return position.x; }

double particle::get_y() { return position.y; }

double particle::get_z() { return position.z; }

void particle::set_position(point position) { this->position = position; }

void particle::set_velocity(point velocity) { this->velocity = velocity; }

void particle::set_x(double x) { this->position.x = x; }

void particle::set_y(double y) { this->position.y = y; }

void particle::set_z(double z) { this->position.z = z; }

void particle::set_u(double u) { this->velocity.x = u; }

void particle::set_v(double v) { this->velocity.y = v; }

void particle::set_w(double w) { this->velocity.z = w; }
}  // namespace mc3d