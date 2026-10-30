#pragma once

#include <functional>
#include <vector>

#include "boundary.h"
#include "utils/utils.hpp"
#include "geometry.h"
#include "poligon.h"

namespace mc3d {
class InnerBoundary {
 private:
  std::vector<std::reference_wrapper<Polygon>> poligon_ptrs;
  double Tw;
  static constexpr double kEpsilon = 0.00001;
  Geometry* geometry_ = nullptr;  // Non-owning; geometry outlives the boundary.
  Point bounds_min_;
  Point bounds_max_;
  bool SegmentMayReachBody(Point start, Point end) const;

 public:
  InnerBoundary();
  InnerBoundary(mc3d::Geometry& geometry, Point cell_center, double L);
  ~InnerBoundary();
  bool AddPolygon(mc3d::Geometry& body, Point cell_center, double L);
  void SetGeometry(Geometry& geometry);
  void SetWallTemperature(double temperature);
  Geometry* GetGeometryPtr();
  bool Empty() const;
  int BoundaryCondition(std::vector<Particle>& cluster_particle, double dt, utils::RandomStream& random);
  double CalcCellVolume(
      Point cell_apex, Point size, utils::RandomStream& random,
      Point* mass_center =
          NULL);  // вычисление отсеченного обёма ячейки, аргументы:
                  // входные параметры:cell_apex - опорный угол ячейки, size -
                  // размеры ячейки, выходные параметры: mass_center - центр мас
                  // ячейки.
};
}  // namespace mc3d
