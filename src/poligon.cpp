// File: Polygon.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of rarefied flows.

#include "poligon.h"

#include <cmath>

#include "exception.h"

namespace mc3d {
Polygon::Polygon() {
  force = Point(0, 0, 0);
  flux = 0;
  this->p1 = Point(0, 0, 0);
  this->p2 = Point(0, 0, 0);
  this->p3 = Point(0, 0, 0);
  this->normal = Point(1, 0, 0);
}

Polygon::~Polygon() {}

Polygon::Polygon(Point p1, Point p2, Point p3, Point nrml) {
  force = Point(0, 0, 0);
  flux = 0;
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  this->normal = nrml;
  S = (p2 - p1).Cross(p3 - p2).Mod();
}

Polygon::Polygon(Point p1, Point p2, Point p3) {
  force = Point(0, 0, 0);
  flux = 0;
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  // this->normal = nrml;
  S = (p2 - p1).Cross(p3 - p2).Mod();
}

void Polygon::Print() {
  p1.Print();
  p2.Print();
  p3.Print();
  normal.Print();
}

void Polygon::Set(Point p1, Point p2, Point p3, Point nrml) {
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  this->normal = nrml;
  S = (p2 - p1).Cross(p3 - p2).Mod();
}

Point Polygon::GetNormal() const { return normal; }

Point Polygon::GetP1() const { return p1; }

Point Polygon::GetP2() const { return p2; }

Point Polygon::GetP3() const { return p3; }

void Polygon::Move(const Point& a) {
  p1 += a;
  p2 += a;
  p3 += a;
}

void Polygon::Scale(double e) {
  p1 *= e;
  p2 *= e;
  p3 *= e;
}

Point Polygon::GetGmt() const { return (p1 + p2 + p3) / 3; }

double Polygon::GetLmax() const {
  double L1 = (p1 - p2).Mod();
  double L2 = (p3 - p2).Mod();
  double L3 = (p1 - p3).Mod();
  double Lmax = L1 > L2 ? L1 : L2;
  Lmax = L3 > Lmax ? L3 : Lmax;
  return Lmax;
}

std::pair<std::unique_ptr<Polygon>, std::unique_ptr<Polygon>> Polygon::Divide()
    const {
  double a1 = (p1 - p2).Mod();
  double a2 = (p2 - p3).Mod();
  double a3 = (p3 - p1).Mod();
  if ((a1 >= a2) && (a1 >= a3)) {
    Point new_pt = (p1 + p2) / 2;
    return {std::make_unique<Polygon>(p1, new_pt, p3, normal),
            std::make_unique<Polygon>(new_pt, p2, p3, normal)};
  }
  if ((a2 >= a1) && (a2 >= a3)) {
    Point new_pt = (p2 + p3) / 2;
    return {std::make_unique<Polygon>(p2, new_pt, p1, normal),
            std::make_unique<Polygon>(new_pt, p3, p1, normal)};
  }
  if ((a3 >= a1) && (a3 >= a2)) {
    Point new_pt = (p3 + p1) / 2;
    return {std::make_unique<Polygon>(p1, p2, new_pt, normal),
            std::make_unique<Polygon>(new_pt, p2, p3, normal)};
  }
  throw unusual_situations::Exception(
      unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR);
  return {std::unique_ptr<Polygon>(), std::unique_ptr<Polygon>()};
}

bool Polygon::Fix() {
  if ((p2 - p1).Cross(p3 - p2) * normal < 0) {
    Point tmp = p3;
    p3 = p2;
    p2 = tmp;
    return false;
  }
  return true;

  /*if((p3 - p2).Cross(p1 - p3) * normal < 0)
  {
          Point tmp = p3;
          p3 = p2;
          p2 = tmp;
  }

  if((p1 - p3).Cross(p2 - p1) * normal < 0)
  {
          Point tmp = p3;
          p3 = p2;
          p2 = tmp;
  }*/
}

double Polygon::DistanceToPoint(Point p) {
  double dist = (p1 - p) * normal;

  Point cpstn = p + normal * dist;

  Point a = p2 - p1;
  Point b = p3 - p2;
  Point c = p1 - p3;

  Point d1 = p1 - cpstn;
  Point d2 = p2 - cpstn;
  Point d3 = p3 - cpstn;

  /*cout << normal * d1.Cross(a) << endl;
  cout << normal * d2.Cross(b) << endl;
  cout << normal * d3.Cross(c) << endl << endl;*/

  if (normal * d1.Cross(a) >= 0. && normal * d2.Cross(b) >= 0. &&
      normal * d3.Cross(c) >= 0.) {
    return std::abs(dist);
  } else {
    return std::min((p1 - p).Mod(), std::min((p3 - p).Mod(), (p3 - p).Mod()));
  }
  return -1.;
}
}  // namespace mc3d
