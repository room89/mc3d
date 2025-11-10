// File: point.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of freemolecular flows.

#pragma once

#include <iostream>

using namespace std;

namespace mc3d {
struct Point {
 public:
  double x, y, z;

 public:
  Point();
  Point(double x, double y, double z);
  ~Point();
  static Point RandomPoint(
      Point size);  // возвращает рандомную точку с коорд.: x от 0 до size.x, y
                    // от 0 до size.y, z от 0 до size.z
  Point
  RandomPoint();  // возвращает рандомную точку с коорд.: x от 0 до this->x,
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
};  // namespace mc3d
