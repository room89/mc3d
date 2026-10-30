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

  Particle() : velocity(0, 0, 0), position(0, 0, 0) {}
  Particle(double x, double y, double z, double u, double v, double w)
      : velocity(u, v, w), position(x, y, z) {}
  Particle(Point position, Point velocity) : velocity(velocity), position(position) {}

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

inline void Particle::Move(double dt) { position += velocity * dt; }

inline double Particle::GetU() const { return velocity.x; }

inline double Particle::GetV() const { return velocity.y; }

inline double Particle::GetW() const { return velocity.z; }

inline double Particle::GetX() const { return position.x; }

inline double Particle::GetY() const { return position.y; }

inline double Particle::GetZ() const { return position.z; }

inline void Particle::SetPosition(Point position) { this->position = position; }

inline void Particle::SetVelocity(Point velocity) { this->velocity = velocity; }

inline void Particle::SetX(double x) { this->position.x = x; }

inline void Particle::SetY(double y) { this->position.y = y; }

inline void Particle::SetZ(double z) { this->position.z = z; }

inline void Particle::SetU(double u) { this->velocity.x = u; }

inline void Particle::SetV(double v) { this->velocity.y = v; }

inline void Particle::SetW(double w) { this->velocity.z = w; }

bool Collision(Particle& a, Particle& b, double& g_max_ch,
               double& frequency_t_ch, double factor, utils::RandomStream& random);
}  // namespace mc3d
