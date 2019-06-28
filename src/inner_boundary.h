#pragma once

#include "boundary.h"
#include "geometry.h"
#include "poligon.h"

namespace mc3d {
class inner_boundary : public boundary {
 private:
  std::deque<poligon*> poligon_ptrs;
  point cell_center;
  double Tw;
  double eps;
  geometry* body;

 public:
  inner_boundary(void);
  inner_boundary(mc3d::geometry* bbody, point cell_center, double L);
  ~inner_boundary(void);
  bool add_poligon(mc3d::geometry* bbody, point cell_center, double L);
  // bool add_poligon(mc3d::geometry *bbody, double L);
  void set_geometry(geometry* bbody);
  geometry* get_geometry_ptr();
  int bondary_condition(deque<particle>* cluster_particle, double dt);
  double calc_cell_volume(
      point cell_apex, point size,
      point* mass_center =
          NULL);  //вычисление отсеченного обёма ячейки, аргументы:
                  //входные параметры:cell_apex - опорный угол ячейки, size -
                  //размеры ячейки, выходные параметры: mass_center - центр мас
                  //ячейки.
};
}  // namespace mc3d