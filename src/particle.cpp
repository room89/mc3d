// File: Particle.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.2.1
// Last modified: 24.12.08.
// Description: Program for calculation of rarefaid flows.

#include "particle.h"

#include <cmath>
#include <numbers>
#include <utils/utils.hpp>

namespace mc3d {
Particle::Particle() {
  // this->r = 0.001;
  // this->u = this->v = this->w = this->x = this->y = this->z = 0;
  position.Set(0, 0, 0);
  velocity.Set(0, 0, 0);
}

Particle::Particle(double x, double y, double z, double u, double v, double w) {
  position.Set(x, y, z);
  velocity.Set(u, v, w);
}

Particle::Particle(Point position, Point velocity) {
  this->position = position;
  this->velocity = velocity;
}

Particle::~Particle() {}

/*void Particle::print()
{
        std::cout<<" x = "<<x<<" y = "<<y<<" z = "<<z<<" u = "<<u<<" v = "<<v<<"
w = "<<w<<std::endl;
}*/

void Particle::Move(double dt) {
  position += velocity * dt;

  /*x += dt * u;
  y += dt * v;
  z += dt * w;*/
}
bool Collision(Particle& a, Particle& b, double& g_max, double& frequency_t,
               double factor) {
  const double Pi = std::numbers::pi;

  /*double u1 = a->u, v1 = a->v, w1 = a->w;
  double u2 = b->u, v2 = b->v, w2 = b->w;
  double gx = u2 - u1, gy = v2 - v1, gz = w2 - w1;
  double g = sqrt(gx * gx + gy * gy + gz * gz);*/

  Point vel_1 = a.GetVelocity();
  Point vel_2 = b.GetVelocity();
  double g = (vel_2 - vel_1).Mod();

  bool rtrn = true;

  if (g > g_max) {
    g_max = g;
    frequency_t = factor / g_max;
    rtrn = false;
  }

  double rr = utils::Random01();

  if (g / g_max > rr) {
    double r1 = utils::Random01();
    double r2 = utils::Random01();

    double g1z = g * std::cos(Pi * r1);
    double g1y = g * std::sin(Pi * r1) * std::sin(2. * Pi * r2);
    double g1x = g * std::sin(Pi * r1) * std::cos(2. * Pi * r2);

    Point g1(g1x, g1y, g1z);

    /*a->u = 0.5 * (u1 + u2) - 0.5 * g1x;
    a->v = 0.5 * (v1 + v2) - 0.5 * g1y;
    a->w = 0.5 * (w1 + w2) - 0.5 * g1z;*/

    a.SetVelocity(0.5 * (vel_1 + vel_2 - g1));

    /*b->u = 0.5 * (u1 + u2) + 0.5 * g1x;
    b->v = 0.5 * (v1 + v2) + 0.5 * g1y;
    b->w = 0.5 * (w1 + w2) + 0.5 * g1z;*/

    b.SetVelocity(0.5 * (vel_1 + vel_2 + g1));
  }

  return rtrn;
}

double Particle::GetU() const { return velocity.x; }

double Particle::GetV() const { return velocity.y; }

double Particle::GetW() const { return velocity.z; }

double Particle::GetX() const { return position.x; }

double Particle::GetY() const { return position.y; }

double Particle::GetZ() const { return position.z; }

void Particle::SetPosition(Point position) { this->position = position; }

void Particle::SetVelocity(Point velocity) { this->velocity = velocity; }

void Particle::SetX(double x) { this->position.x = x; }

void Particle::SetY(double y) { this->position.y = y; }

void Particle::SetZ(double z) { this->position.z = z; }

void Particle::SetU(double u) { this->velocity.x = u; }

void Particle::SetV(double v) { this->velocity.y = v; }

void Particle::SetW(double w) { this->velocity.z = w; }
}  // namespace mc3d