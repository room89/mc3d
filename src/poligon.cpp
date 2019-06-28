// File: poligon.cpp
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Last modified: 12.10.10.
// Description: Program for calculation of rarefied flows.

#include <cmath>

#include "poligon.h"

namespace mc3d {
poligon::poligon(void) {
  this->p1 = point(0, 0, 0);
  this->p2 = point(0, 0, 0);
  this->p3 = point(0, 0, 0);
  this->normal = point(1, 0, 0);
}

poligon::~poligon(void) {}

poligon::poligon(point p1, point p2, point p3, point nrml) {
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  this->normal = nrml;
  S = (p2 - p1).vec_mult(p3 - p2).mod();
}

poligon::poligon(point p1, point p2, point p3) {
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  // this->normal = nrml;
  S = (p2 - p1).vec_mult(p3 - p2).mod();
}

void poligon::print() {
  p1.print();
  p2.print();
  p3.print();
  normal.print();
}

void poligon::set(point p1, point p2, point p3, point nrml) {
  this->p1 = p1;
  this->p2 = p2;
  this->p3 = p3;
  this->normal = nrml;
  S = (p2 - p1).vec_mult(p3 - p2).mod();
}

point poligon::get_normal() { return normal; }

point poligon::get_p1() const { return p1; }

point poligon::get_p2() const { return p2; }

point poligon::get_p3() { return p3; }

void poligon::move(point a) {
  p1 += a;
  p2 += a;
  p3 += a;
}

void poligon::scaling(double e) {
  p1 *= e;
  p2 *= e;
  p3 *= e;
}

point poligon::get_gmt() const { return (p1 + p2 + p3) / 3; }

poligon* poligon::get_ptr() { return this; }

double poligon::get_Lmax() {
  double L1 = (p1 - p2).mod();
  double L2 = (p3 - p2).mod();
  double L3 = (p1 - p3).mod();
  double Lmax = L1 > L2 ? L1 : L2;
  Lmax = L3 > Lmax ? L3 : Lmax;
  return Lmax;
}

pair<poligon*, poligon*> poligon::divide() {
  double a1 = (p1 - p2).mod();
  double a2 = (p2 - p3).mod();
  double a3 = (p3 - p1).mod();
  if ((a1 >= a2) && (a1 >= a3)) {
    point new_pt = (p1 + p2) / 2;
    pair<poligon*, poligon*> new_poligon_pair;
    new_poligon_pair.first = new poligon(p1, new_pt, p3, normal);
    new_poligon_pair.second = new poligon(new_pt, p2, p3, normal);
    return new_poligon_pair;
  }
  if ((a2 >= a1) && (a2 >= a3)) {
    point new_pt = (p2 + p3) / 2;
    pair<poligon*, poligon*> new_poligon_pair;
    new_poligon_pair.first = new poligon(p2, new_pt, p1, normal);
    new_poligon_pair.second = new poligon(new_pt, p3, p1, normal);
    return new_poligon_pair;
  }
  if ((a3 >= a1) && (a3 >= a2)) {
    point new_pt = (p3 + p1) / 2;
    pair<poligon*, poligon*> new_poligon_pair;
    new_poligon_pair.first = new poligon(p1, p2, new_pt, normal);
    new_poligon_pair.second = new poligon(new_pt, p2, p3, normal);
    return new_poligon_pair;
  }
  throw unusual_situations::exception(
      unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR, this);
  poligon* null_poligon_ptr = NULL;
  pair<poligon*, poligon*> error_poligon_pair(null_poligon_ptr,
                                              null_poligon_ptr);
  return error_poligon_pair;
}

bool poligon::fix() {
  if ((p2 - p1).vec_mult(p3 - p2) * normal < 0) {
    point tmp = p3;
    p3 = p2;
    p2 = tmp;
    return false;
  }
  return true;

  /*if((p3 - p2).vec_mult(p1 - p3) * normal < 0)
  {
          point tmp = p3;
          p3 = p2;
          p2 = tmp;
  }

  if((p1 - p3).vec_mult(p2 - p1) * normal < 0)
  {
          point tmp = p3;
          p3 = p2;
          p2 = tmp;
  }*/
}

double poligon::dist_to_point(point p) {
  double dist = (p1 - p) * normal;

  point cpstn = p + normal * dist;

  point a = p2 - p1;
  point b = p3 - p2;
  point c = p1 - p3;

  point d1 = p1 - cpstn;
  point d2 = p2 - cpstn;
  point d3 = p3 - cpstn;

  /*cout << normal * d1.vec_mult(a) << endl;
  cout << normal * d2.vec_mult(b) << endl;
  cout << normal * d3.vec_mult(c) << endl << endl;*/

  if (normal * d1.vec_mult(a) >= 0. && normal * d2.vec_mult(b) >= 0. &&
      normal * d3.vec_mult(c) >= 0.) {
    return std::abs(dist);
  } else {
    return std::min((p1 - p).mod(), std::min((p3 - p).mod(), (p3 - p).mod()));
  }
  return -1.;
}
}  // namespace mc3d
