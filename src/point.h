// File: point.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#pragma once

#include <iostream>
#include <cmath>
#include "utils/utils.hpp"

using namespace std;

namespace mc3d {
struct Point {
 public:
  double x, y, z;

 public:
  Point();
  Point(double x, double y, double z);
  static Point RandomPoint(
      Point size, utils::RandomStream& random);  // возвращает рандомную точку с коорд.: x от 0 до size.x, y
                    // от 0 до size.y, z от 0 до size.z
  Point
  RandomPoint(utils::RandomStream& random);  // возвращает рандомную точку с коорд.: x от 0 до this->x,
                  //  y от 0 до this->y, z от 0 до this->z
  void Print();
  double Mod() const;
  double Volume() const;
  Point& Normalize();
  Point Set(double x, double y, double z);
  double GetX();
  Point Get();
  Point operator+=(Point b);
  Point operator-=(Point b);
  Point operator/=(double b);
  Point operator*=(double b);
  // Point& operator=(const Point &a);
  friend Point operator+(const Point a, const Point b);
  friend double operator*(const Point a, const Point b);
  friend Point operator*(const Point a, const double b);
  friend Point operator*(const double a, const Point b);
  friend Point operator/(const Point a, const double b);
  friend Point operator-(const Point a, const Point b);
  friend Point operator!(const Point a);
  friend ostream& operator<<(ostream&, Point);
  Point Cross(Point p1);
};

// Keep default/value initialization unchanged: coordinates are not initialized.
inline Point::Point() {}

inline Point::Point(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

inline Point Point::Set(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
  return Point(x, y, z);
}

inline Point operator+(const Point a, const Point b) {
  return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

inline double operator*(const Point a, const Point b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline Point operator*(Point a, double b) { return Point(a.x * b, a.y * b, a.z * b); }

inline Point operator*(double a, Point b) { return Point(b.x * a, b.y * a, b.z * a); }

inline Point operator/(Point a, double b) { return a.Set(a.x / b, a.y / b, a.z / b); }

inline Point operator-(const Point a, const Point b) {
  return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

inline Point operator!(Point a) {
  return Point(a.x < 0.0000001 ? 1 : 0, a.y < 0.0000001 ? 1 : 0,
               a.z < 0.0000001 ? 1 : 0);
}

inline double Point::Mod() const { return std::sqrt(x * x + y * y + z * z); }

inline Point Point::Cross(Point a) {
  return Point(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}

inline Point& Point::Normalize() {
  double m = Mod();
  Point a(x / m, y / m, z / m);
  *this = a;
  return *this;
}

inline double Point::GetX() { return x; }

inline Point Point::operator+=(Point a) {
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}

inline Point Point::operator-=(Point a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}

inline Point Point::operator/=(double a) {
  x /= a;
  y /= a;
  z /= a;
  return *this;
}

inline Point Point::operator*=(double a) {
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

inline Point Point::Get() {
  Point a;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

inline double Point::Volume() const { return x * y * z; }

}  // namespace mc3d
