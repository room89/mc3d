#pragma once

#include "boundary.h"
#include "geometry.h"
#include "poligon.h"

namespace mc3d {
class InnerBoundary : public Boundary {
 private:
  std::deque<Polygon*> poligon_ptrs;
  Point cell_center;
  double Tw;
  double eps;
  Geometry* body;

 public:
  InnerBoundary();
  InnerBoundary(mc3d::Geometry* bbody, Point cell_center, double L);
  ~InnerBoundary();
  bool AddPolygon(mc3d::Geometry* bbody, Point cell_center, double L);
  // bool add_poligon(mc3d::Geometry *bbody, double L);
  void SetGeometry(Geometry* bbody);
  Geometry* GetGeometryPtr();
  bool Empty() const;
  int BoundaryCondition(deque<Particle>& cluster_particle, double dt);
  double CalcCellVolume(
      Point cell_apex, Point size,
      Point* mass_center =
          NULL);  // вычисление отсеченного обёма ячейки, аргументы:
                  // входные параметры:cell_apex - опорный угол ячейки, size -
                  // размеры ячейки, выходные параметры: mass_center - центр мас
                  // ячейки.
};
}  // namespace mc3d
