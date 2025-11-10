// File: point.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#include "point.h"

#include <cmath>
#include <iostream>
#include <utils/utils.hpp>

namespace mc3d {
Point::Point() {}

Point::~Point() {}

Point::Point(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

void Point::Print() {
  std::cout << "(" << x << "," << y << "," << this->z << ")" << std::endl;
}

Point Point::Set(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
  return Point(x, y, z);
}

Point operator+(const Point a, const Point b) {
  return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

double operator*(const Point a, const Point b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

Point operator*(Point a, double b) { return Point(a.x * b, a.y * b, a.z * b); }

Point operator*(double a, Point b) { return Point(b.x * a, b.y * a, b.z * a); }

Point operator/(Point a, double b) { return a.Set(a.x / b, a.y / b, a.z / b); }

Point operator-(const Point a, const Point b) {
  return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

Point operator!(Point a) {
  return Point(a.x < 0.0000001 ? 1 : 0, a.y < 0.0000001 ? 1 : 0,
               a.z < 0.0000001 ? 1 : 0);
}

std::ostream& operator<<(std::ostream& o, const Point c) {
  o << c.x << ", " << c.y << ", " << c.z;
  return o;
}

double Point::Mod() const { return std::sqrt(x * x + y * y + z * z); }

Point Point::Cross(Point a) {
  return Point(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}

Point& Point::Normalize() {
  double m = Mod();
  Point a(x / m, y / m, z / m);
  *this = a;
  return *this;
}

double Point::GetX() { return x; }

Point Point::operator+=(Point a) {
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}

Point Point::operator-=(Point a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}

Point Point::operator/=(double a) {
  x /= a;
  y /= a;
  z /= a;
  return *this;
}

Point Point::operator*=(double a) {
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

/*Point& Point::operator=(const Point &a)
{
        x = a.x;
        y = a.y;
        z = a.z;
        return *this;
}*/

Point Point::Get() {
  Point a;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

Point Point::RandomPoint(
    Point size)  // возврощает рандомную точку с коорд.: x от 0 до size.x, y от
                 // 0 до size.y, z от 0 до size.z
{
  double rnx = utils::RandomDouble(0.0, 1.0);
  double rny = utils::RandomDouble(0.0, 1.0);
  double rnz = utils::RandomDouble(0.0, 1.0);

  return Point(size.x * rnx, size.y * rny, size.z * rnz);
}

Point Point::RandomPoint()  // возврощает рандомную точку с коорд.: x от 0 до
                            //  this->x, y от 0 до this->y, z от 0 до this->z
{
  double rnx = utils::RandomDouble(0.0, 1.0);
  double rny = utils::RandomDouble(0.0, 1.0);
  double rnz = utils::RandomDouble(0.0, 1.0);

  return Point(x * rnx, y * rny, z * rnz);
}

double Point::Volume() const { return x * y * z; }
}  // namespace mc3d
