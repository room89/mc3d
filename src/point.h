// File: point.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#pragma once

#include <iostream>

using namespace std;

namespace mc3d {
struct point {
 public:
  double x, y, z;

 public:
  point(void);
  point(double x, double y, double z);
  ~point(void);
  static point rand_point(
      point size);  //возврощает рандомную точку с коорд.: x от 0 до size.x, y
                    //от 0 до size.y, z от 0 до size.z
  point rand_point();  //возврощает рандомную точку с коорд.: x от 0 до this->x,
                       //y от 0 до this->y, z от 0 до this->z
  void print();
  double mod() const;
  double volume() const;
  point nrmlz();
  point set(double x, double y, double z);
  double get_x();
  point get();
  point operator+=(point b);
  point operator-=(point b);
  point operator/=(double b);
  point operator*=(double b);
  // point& operator=(const point &a);
  friend point operator+(const point a, const point b);
  friend double operator*(const point a, const point b);
  friend point operator*(const point a, const double b);
  friend point operator*(const double a, const point b);
  friend point operator/(const point a, const double b);
  friend point operator-(const point a, const point b);
  friend point operator!(const point a);
  friend ostream& operator<<(ostream&, point);
  point vec_mult(point p1);
};
};  // namespace mc3d