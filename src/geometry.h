// File: Geometry.h
// Program: MC3D
// Author: Khokhlov "AAsad" Ivan
// Version: 0.4.0
// Last modified: 11.10.10.
// Description: Program for calculation of rarefied flows.

#pragma once

#include <cstddef>
#include <deque>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "boundary.h"
#include "exception.h"
#include "poligon.h"
// #include "inner_boundary.h"

namespace mc3d {
namespace {
// const double Pi = 3.141592654;
}

class Geometry {
 protected:
  std::deque<std::unique_ptr<Polygon>> poligons;
  std::mutex force_mutex_;
  // Subclasses editing polygon vertices/topology must invalidate before edits.
  void InvalidateQueries() {
    query_bounds_.reset();
    query_tolerance_.reset();
    query_triangles_.clear();
  }

 private:
  struct QueryTriangle {
    Point origin;
    Point e1;
    Point e2;
    double determinant_scale;
  };
  static QueryTriangle PrepareTriangle(const Polygon& polygon);
  std::optional<double> RayTriangle(std::size_t index, Point start,
                                    Point direction, double direction_length) const;
  static std::optional<double> IntersectTriangle(const QueryTriangle& triangle,
      Point start, Point direction, double direction_length);
  std::vector<QueryTriangle> query_triangles_;
  std::optional<std::pair<Point, Point>> query_bounds_;
  std::optional<double> query_tolerance_;

 public:
  struct SurfaceHit {
    Polygon* polygon;
    double fraction;
  };
  Geometry();
  Geometry(const char* file_name);
  ~Geometry();
  void WriteGeometryFile(const char* file_name, const char* boby_name);
  void Fragment(double Lmax);
  // void CreateCone(Point apex1, Point apex2, double H);
  void CreateWedge(double x = -.25, double width = .5, double length = .5,
                   double alpha = 3.141592654 / 18);
  void CreatePyramid(double x = -.25, double width = .5, double length = .5,
                     double H = .16);
  void CreateCube(double x = -.25, double width = .5, double length = .5,
                  double H = .5);
  void CreateCylinder(double x = -.25, double radius = .25, double length = .5,
                      int segments = 32);
  Point MassCenter();
  void Move(Point a);
  void Scale(double e);
  std::size_t PolygonCount() const;
  const Polygon& GetPolygon(std::size_t index) const;
  bool IsInnerPoint(Point test_point) const;
  std::optional<SurfaceHit> FirstIntersection(Point start, Point displacement);
  std::optional<Point> ExteriorPoint(Point interior) const;
  // Prepare on the setup thread. Geometry must not be edited during queries.
  void PrepareQueries();
  std::pair<Point, Point> Bounds() const;
  double SurfaceTolerance() const;
  void AccumulateForce(Polygon& polygon, Point impulse);
  pair<Point, Point> Size();
  double DistanceToPoint(Point p);
  int FixPolygons();

  void ReverseNormals();

  friend class InnerBoundary;
};
}  // namespace mc3d
