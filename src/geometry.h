// File: geometry.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.4.0
// Last modified: 11.10.10.
// Description: Program for calculation of rarefied flows.

#pragma once

#include <deque>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include "boundary.h"
#include "exception.h"
#include "poligon.h"
//#include "inner_boundary.h"

namespace mc3d {
namespace {
// const double Pi = 3.141592654;
}

class geometry {
 protected:
  deque<poligon*> poligons;

 public:
  geometry(void);
  geometry(const char* file_name);
  ~geometry(void);
  void write_geometry_file(const char* file_name, const char* boby_name);
  void fragmentation(double Lmax);
  // void create_cone(point apex1, point apex2, double H);
  void create_wedge(double x = -.25, double width = .5, double length = .5,
                    double alpha = 3.141592654 / 18);
  void create_pyramid(double x = -.25, double width = .5, double length = .5,
                      double H = .16);
  void create_cube(double x = -.25, double width = .5, double length = .5,
                   double H = .5);
  point mass_center();
  void move(point a);
  void scaling(double e);
  bool is_inner_point(point test_point);
  pair<point, point> size();
  double dist_to_point(point p);
  int fix_poligons();

  void rev_nrml();

  friend class inner_boundary;
};
}  // namespace mc3d
