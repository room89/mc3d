#include "inner_boundary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utils/utils.hpp>

namespace mc3d {
const double Pi = 3.1415926535;

InnerBoundary::InnerBoundary() {
  Tw = 1;
  eps = 0.00001;
}

InnerBoundary::InnerBoundary(mc3d::Geometry& geometry, Point cell_center,
                             double L)
    : InnerBoundary() {
  SetGeometry(geometry);
  AddPolygon(geometry, cell_center, L);
}

InnerBoundary::~InnerBoundary() {}

bool InnerBoundary::AddPolygon(Geometry& body, Point cell_center, double L) {
  this->cell_center = cell_center;

  poligon_ptrs.clear();

  auto poligon_iter = body.poligons.begin();

  while (poligon_iter != body.poligons.end()) {
    Polygon& polygon = *(*poligon_iter);
    polygon.Fix();
    const auto& distance_to_gmt = (polygon.GetGmt() - cell_center).Mod();
    const auto& distance_to_p1 = (polygon.GetP1() - cell_center).Mod();
    const auto& distance_to_p2 = (polygon.GetP2() - cell_center).Mod();
    const auto& distance_to_p3 = (polygon.GetP3() - cell_center).Mod();
    if (distance_to_gmt < L    // проверка растояния до ГМТ полигона
        || distance_to_p1 < L  // проверка растояния до первого угла полигона
        || distance_to_p2 < L  // проверка растояния до второго угла полигона
        || distance_to_p3 < L  // проверка растояния до третьего угла полигона
    ) {
      poligon_ptrs.emplace_back(polygon);
    } else {
      double dist = (cell_center - polygon.GetGmt()) * polygon.normal;
      if ((dist < L) && (dist > -L)) {
        Point point_on_poligon = cell_center + polygon.normal * dist;

        polygon.Fix();

        Point a = polygon.GetP2() - polygon.GetP1();
        Point b = polygon.GetP3() - polygon.GetP2();
        Point c = polygon.GetP1() - polygon.GetP3();

        Point d1 = polygon.GetP1() - point_on_poligon;
        Point d2 = polygon.GetP2() - point_on_poligon;
        Point d3 = polygon.GetP3() - point_on_poligon;

        // Более строгая проверка принадлежности точки полигону
        if (polygon.GetNormal() * d1.Cross(a) > -eps * 10.0 &&
            polygon.GetNormal() * d2.Cross(b) > -eps * 10.0 &&
            polygon.GetNormal() * d3.Cross(c) > -eps * 10.0) {
          poligon_ptrs.emplace_back(polygon);
        }
      }
    }

    poligon_iter++;
  }

  return !poligon_ptrs.empty();
}

int InnerBoundary::BoundaryCondition(std::vector<Particle>& cluster_particle,
                                     double dt) {
  Geometry* geometry_ptr = GetGeometryPtr();
  if (!geometry_ptr) return 0;
  const double surface_offset = 8 * geometry_ptr->SurfaceTolerance();

  for (auto& particle : cluster_particle) {
    if (!SegmentMayReachBody(particle.position,
                             particle.position + particle.velocity * dt)) {
      particle.position += particle.velocity * dt;
      continue;
    }
    if (geometry_ptr->IsInnerPoint(particle.position)) {
      if (auto exterior = geometry_ptr->ExteriorPoint(particle.position)) {
        particle.position = *exterior;
      } else {
        throw std::runtime_error("Cannot recover particle from body interior");
      }
    }

    double particle_dt = dt;
    for (int collisions = 0; particle_dt > 0 && collisions < 32; ++collisions) {
      if (!SegmentMayReachBody(
              particle.position,
              particle.position + particle.velocity * particle_dt)) {
        particle.position += particle.velocity * particle_dt;
        particle_dt = 0;
        break;
      }
      const auto hit = geometry_ptr->FirstIntersection(
          particle.position, particle.velocity * particle_dt);
      if (!hit) {
        particle.position += particle.velocity * particle_dt;
        particle_dt = 0;
        break;
      }

      Polygon* collision_polygon = hit->polygon;
      const double dtt = particle_dt * hit->fraction;
      Point edge1 = collision_polygon->GetP2() - collision_polygon->GetP1();
      Point edge2 = collision_polygon->GetP3() - collision_polygon->GetP1();
      Point unit_normal = edge1.Cross(edge2);
      const double normal_length = unit_normal.Mod();
      if (!std::isfinite(normal_length) || normal_length == 0) {
        throw std::runtime_error("Invalid triangle at particle collision");
      }
      unit_normal /= normal_length;
      // The mesh normal may be reversed; an entering particle must leave the
      // collision along the side from which it approached the triangle.
      if (unit_normal * particle.velocity > 0) unit_normal *= -1;
      particle.position +=
          particle.velocity * dtt + unit_normal * surface_offset;

      const Point incoming_velocity = particle.velocity;

      double r1 = utils::Random01();
      double r2 = utils::Random01();
      double r3 = utils::Random01();

      if (r1 <= 0.) r1 = 0.00001;
      if (r2 <= 0.) r2 = 0.00001;

      Point vel = particle.velocity;

      double sp = vel * unit_normal;
      double r = sqrt(2. * Tw * fabs(log(r1)));

      Point vn;
      vn = unit_normal * r;

      Point vni;
      vni = unit_normal * sp;

      Point vt;
      vt = vel - vni;
      double lvt = vt.Mod();
      Point tangent1;
      if (lvt < eps) {
        if (std::abs(unit_normal.x) < 0.9) {
          tangent1 = unit_normal.Cross(Point(1.0, 0.0, 0.0));
        } else {
          tangent1 = unit_normal.Cross(Point(0.0, 1.0, 0.0));
        }
        double tangent1_length = tangent1.Mod();
        if (tangent1_length < eps) {
          tangent1 = unit_normal.Cross(Point(0.0, 0.0, 1.0));
          tangent1_length = tangent1.Mod();
          if (tangent1_length < eps) {
            tangent1 = Point(0.0, 1.0, 0.0);
            tangent1_length = tangent1.Mod();
          }
        }
        tangent1 /= tangent1_length;
      } else {
        vt /= lvt;
        tangent1 = vt;
      }

      r = sqrt(2. * Tw * fabs(log(r2)));
      double teta = 2. * Pi * r3;
      double vt1m = r * cos(teta);
      double vt2m = r * sin(teta);

      Point vt1 = tangent1;
      Point vt2;

      vt1 *= vt1m;

      vt2 = unit_normal.Cross(tangent1);
      double tangent2_length = vt2.Mod();
      if (tangent2_length < eps) {
        if (std::abs(unit_normal.z) < 0.9) {
          vt2 = unit_normal.Cross(Point(0.0, 0.0, 1.0));
        } else {
          vt2 = unit_normal.Cross(Point(0.0, 1.0, 0.0));
        }
        tangent2_length = vt2.Mod();
        if (tangent2_length < eps) {
          vt2 = Point(1.0, 0.0, 0.0);
          tangent2_length = vt2.Mod();
        }
      }
      vt2 /= tangent2_length;
      vt2 *= vt2m;

      particle.velocity = vt1 + vt2 + vn;
      double outgoing_component = particle.velocity * unit_normal;
      if (outgoing_component <= 0.0) {
        particle.velocity -= unit_normal * (2.0 * outgoing_component);
      }

      geometry_ptr->AccumulateForce(*collision_polygon,
                                    incoming_velocity - particle.velocity);

      if (geometry_ptr->IsInnerPoint(particle.position)) {
        if (auto exterior = geometry_ptr->ExteriorPoint(particle.position)) {
          particle.position = *exterior;
        } else {
          throw std::runtime_error(
              "Cannot recover particle from body interior");
        }
      }
      particle_dt -= dtt;
    }
    if (particle_dt > 0) {
      throw std::runtime_error(
          "Body collision limit reached with unprocessed step time; reduce Cu");
    }
    if (geometry_ptr->IsInnerPoint(particle.position)) {
      if (auto exterior = geometry_ptr->ExteriorPoint(particle.position)) {
        particle.position = *exterior;
      } else {
        throw std::runtime_error("Cannot recover particle from body interior");
      }
    }
  }

  return 0;
}

double InnerBoundary::CalcCellVolume(
    Point cell_apex, Point size,
    Point* mass_center_out)  // вычисление отсеченного обёма ячейки, аргументы:
                             // входные параметры:cell_apex - опорный угол
                             // ячейки, size - размеры ячейки, выходные
                             // параметры: mass_center - центр мас ячейки.
{
  double volume = 0;
  double tot_vol = size.x * size.y * size.z;
  double dtx = size.x / 1;
  double dty = size.y / 1;
  double dtz = size.z / 1;

  Point mass_center(0, 0, 0);

  if (poligon_ptrs.size() == 0) {
    // Заглушка!!!
    return tot_vol;
  }

  size_t N = 100000;
  size_t outer_N = N;

  std::vector<Point> rand_points;
  rand_points.reserve(N);

  for (size_t i = 0; i < N; i++) {
    rand_points.emplace_back(cell_apex + size.RandomPoint());
  }

  auto poligon_col_iterator_x = poligon_ptrs.end();
  auto poligon_col_iterator_y = poligon_ptrs.end();
  auto poligon_col_iterator_z = poligon_ptrs.end();

  Point collision_pstn_x;
  Point collision_pstn_y;
  Point collision_pstn_z;

  for (size_t i = 0; i < N; i++) {
    double dttx = size.x / 1;
    double dtty = size.y / 1;
    double dttz = size.z / 1;

    auto poligon_iterator = poligon_ptrs.begin();

    while (poligon_iterator !=
           poligon_ptrs.end())  // ищем полигон с которым будет соударятся
                                // виртуальная частица
    {
      Point poligon_gmt = poligon_iterator->get().GetGmt();

      Point a =
          poligon_iterator->get().GetP2() - poligon_iterator->get().GetP1();
      Point b =
          poligon_iterator->get().GetP3() - poligon_iterator->get().GetP2();
      Point c =
          poligon_iterator->get().GetP1() - poligon_iterator->get().GetP3();

      double Ax = Point(1, 0, 0) * poligon_iterator->get().GetNormal();
      double tcx = (poligon_iterator->get().GetP1() - rand_points[i]) *
                   poligon_iterator->get().GetNormal() / Ax;

      if (tcx < dttx && tcx > -dttx) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcx;

        Point d1 = poligon_iterator->get().GetP1() - collision_pstn;
        Point d2 = poligon_iterator->get().GetP2() - collision_pstn;
        Point d3 = poligon_iterator->get().GetP3() - collision_pstn;

        if (poligon_iterator->get().GetNormal() * d1.Cross(a) > 0. &&
            poligon_iterator->get().GetNormal() * d2.Cross(b) > 0. &&
            poligon_iterator->get().GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_x = poligon_iterator;
          collision_pstn_x = collision_pstn;
          dttx = tcx;
        }
      }

      double Ay = Point(0, 1, 0) * poligon_iterator->get().GetNormal();
      double tcy = (poligon_iterator->get().GetP1() - rand_points[i]) *
                   poligon_iterator->get().GetNormal() / Ay;

      if (tcy < dtty && tcy > -dtty) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcy;

        Point d1 = poligon_iterator->get().GetP1() - collision_pstn;
        Point d2 = poligon_iterator->get().GetP2() - collision_pstn;
        Point d3 = poligon_iterator->get().GetP3() - collision_pstn;

        if (poligon_iterator->get().GetNormal() * d1.Cross(a) > 0. &&
            poligon_iterator->get().GetNormal() * d2.Cross(b) > 0. &&
            poligon_iterator->get().GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_y = poligon_iterator;
          collision_pstn_y = collision_pstn;
          dtty = tcy;
        }
      }

      double Az = Point(0, 0, 1) * poligon_iterator->get().GetNormal();
      double tcz = (poligon_iterator->get().GetP1() - rand_points[i]) *
                   poligon_iterator->get().GetNormal() / Az;

      if (tcz < dttz && tcz > -dttz) {
        Point collision_pstn = rand_points[i] + Point(1, 0, 0) * tcz;

        Point d1 = poligon_iterator->get().GetP1() - collision_pstn;
        Point d2 = poligon_iterator->get().GetP2() - collision_pstn;
        Point d3 = poligon_iterator->get().GetP3() - collision_pstn;

        if (poligon_iterator->get().GetNormal() * d1.Cross(a) > 0. &&
            poligon_iterator->get().GetNormal() * d2.Cross(b) > 0. &&
            poligon_iterator->get().GetNormal() * d3.Cross(c) > 0.) {
          poligon_col_iterator_z = poligon_iterator;
          collision_pstn_z = collision_pstn;
          dttz = tcz;
        }
      }
      poligon_iterator++;
    }  // while(poligon_iterator != poligon_ptrs.end())

    if (poligon_col_iterator_x != poligon_ptrs.end()) {
      if (poligon_col_iterator_x->get().GetNormal() *
              (collision_pstn_x - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_y != poligon_ptrs.end()) {
      if (poligon_col_iterator_y->get().GetNormal() *
              (collision_pstn_y - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    if (poligon_col_iterator_z != poligon_ptrs.end()) {
      if (poligon_col_iterator_z->get().GetNormal() *
              (collision_pstn_z - rand_points[i]) <=
          0) {
        outer_N--;
        break;
      }
    }

    mass_center += rand_points[i];
  }

  volume = tot_vol * (double(outer_N) / double(N));
  if (mass_center_out != NULL) {
    mass_center /= outer_N;
    *mass_center_out = mass_center;
  }

  if (outer_N == N) {
    return -1;
  }

  return volume;
}

void InnerBoundary::SetGeometry(Geometry& geometry) {
  geometry_ = geometry;
  const auto bounds = geometry.Bounds();
  bounds_min_ = bounds.first;
  bounds_max_ = bounds.second;
}

bool InnerBoundary::SegmentMayReachBody(Point start, Point end) const {
  return std::max(start.x, end.x) >= bounds_min_.x &&
         std::min(start.x, end.x) <= bounds_max_.x &&
         std::max(start.y, end.y) >= bounds_min_.y &&
         std::min(start.y, end.y) <= bounds_max_.y &&
         std::max(start.z, end.z) >= bounds_min_.z &&
         std::min(start.z, end.z) <= bounds_max_.z;
}

Geometry* InnerBoundary::GetGeometryPtr() {
  return geometry_ ? &geometry_->get() : nullptr;
}

bool InnerBoundary::Empty() const { return poligon_ptrs.empty(); }
}  // namespace mc3d
