// File: point.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#include "point.h"
#include <complex>
#include <fstream>
#include <iostream>

namespace mc3d {
point::point(void) {}

point::~point(void) {}

point::point(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

void point::print() {
  std::cout << "(" << x << "," << y << "," << this->z << ")" << std::endl;
}

point point::set(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
  return point(x, y, z);
}

point operator+(const point a, const point b) {
  return point(a.x + b.x, a.y + b.y, a.z + b.z);
}

double operator*(const point a, const point b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

point operator*(point a, double b) { return point(a.x * b, a.y * b, a.z * b); }

point operator*(double a, point b) { return point(b.x * a, b.y * a, b.z * a); }

point operator/(point a, double b) { return a.set(a.x / b, a.y / b, a.z / b); }

point operator-(const point a, const point b) {
  return point(a.x - b.x, a.y - b.y, a.z - b.z);
}

point operator!(point a) {
  return point(a.x < 0.0000001 ? 1 : 0, a.y < 0.0000001 ? 1 : 0,
               a.z < 0.0000001 ? 1 : 0);
}

std::ostream& operator<<(std::ostream& o, const point c) {
  o << c.x << "\t" << c.y << "\t" << c.z << "\t";
  return o;
}

double point::mod() const { return std::sqrt(x * x + y * y + z * z); }

point point::vec_mult(point a) {
  return point(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x);
}

point point::nrmlz() {
  double m = mod();
  point a(x / m, y / m, z / m);
  *this = a;
  return a;
}

double point::get_x() { return x; }

point point::operator+=(point a) {
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}

point point::operator-=(point a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}

point point::operator/=(double a) {
  x /= a;
  y /= a;
  z /= a;
  return *this;
}

point point::operator*=(double a) {
  x *= a;
  y *= a;
  z *= a;
  return *this;
}

/*point& point::operator=(const point &a)
{
        x = a.x;
        y = a.y;
        z = a.z;
        return *this;
}*/

point point::get() {
  point a;
  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

point point::rand_point(
    point size)  //возврощает рандомную точку с коорд.: x от 0 до size.x, y от 0
                 //до size.y, z от 0 до size.z
{
  double rmt = 1 / double(RAND_MAX);

  double rnx = std::rand() * rmt;
  double rny = std::rand() * rmt;
  double rnz = std::rand() * rmt;

  return point(size.x * rnx, size.y * rny, size.z * rnz);
}

point point::rand_point()  //возврощает рандомную точку с коорд.: x от 0 до
                           // this->x, y от 0 до this->y, z от 0 до this->z
{
  double rmt = 1 / double(RAND_MAX);

  double rnx = std::rand() * rmt;
  double rny = std::rand() * rmt;
  double rnz = std::rand() * rmt;

  return point(x * rnx, y * rny, z * rnz);
}

double point::volume() const { return x * y * z; }
}  // namespace mc3d