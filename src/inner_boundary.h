#pragma once

#include <functional>
#include <optional>
#include <vector>

#include "boundary.h"
#include "geometry.h"
#include "poligon.h"

namespace mc3d {
class InnerBoundary : public Boundary {
 private:
  std::vector<std::reference_wrapper<Polygon>> poligon_ptrs;
  Point cell_center;
  double Tw;
  double eps;
  std::optional<std::reference_wrapper<Geometry>> geometry_;
  Point bounds_min_;
  Point bounds_max_;
  bool SegmentMayReachBody(Point start, Point end) const;

 public:
  InnerBoundary();
  InnerBoundary(mc3d::Geometry& geometry, Point cell_center, double L);
  ~InnerBoundary();
  bool AddPolygon(mc3d::Geometry& body, Point cell_center, double L);
  void SetGeometry(Geometry& geometry);
  Geometry* GetGeometryPtr();
  bool Empty() const;
  int BoundaryCondition(std::vector<Particle>& cluster_particle, double dt);
  double CalcCellVolume(
      Point cell_apex, Point size,
      Point* mass_center =
          NULL);  // вычисление отсеченного обёма ячейки, аргументы:
                  // входные параметры:cell_apex - опорный угол ячейки, size -
                  // размеры ячейки, выходные параметры: mass_center - центр мас
                  // ячейки.
};
}  // namespace mc3d
