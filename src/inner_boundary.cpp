#include "inner_boundary.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utils/utils.hpp>

namespace mc3d {
const double Pi = 3.1415926535;

InnerBoundary::InnerBoundary() {
  Tw = 1;
}

InnerBoundary::InnerBoundary(mc3d::Geometry& geometry, Point cell_center,
                             double L)
    : InnerBoundary() {
  SetGeometry(geometry);
  AddPolygon(geometry, cell_center, L);
}

InnerBoundary::~InnerBoundary() {}

void InnerBoundary::SetWallTemperature(double temperature) {
  if (!(temperature > 0) || !std::isfinite(temperature))
    throw std::invalid_argument("Wall temperature must be finite and positive");
  Tw = temperature;
}

bool InnerBoundary::AddPolygon(Geometry& body, Point cell_center, double L) {

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
        if (polygon.GetNormal() * d1.Cross(a) > -kEpsilon * 10.0 &&
            polygon.GetNormal() * d2.Cross(b) > -kEpsilon * 10.0 &&
            polygon.GetNormal() * d3.Cross(c) > -kEpsilon * 10.0) {
          poligon_ptrs.emplace_back(polygon);
        }
      }
    }

    poligon_iter++;
  }

  return !poligon_ptrs.empty();
}

int InnerBoundary::BoundaryCondition(std::vector<Particle>& cluster_particle,
                                     double dt, utils::RandomStream& random) {
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
      const double wall_temperature = geometry_ptr->SurfaceHeatingEnabled()
          ? collision_polygon->wall_temperature : Tw;

      double r1 = random.Uniform();
      double r2 = random.Uniform();
      double r3 = random.Uniform();

      if (r1 <= 0.) r1 = 0.00001;
      if (r2 <= 0.) r2 = 0.00001;

      Point vel = particle.velocity;

      double sp = vel * unit_normal;
      double r = sqrt(2. * wall_temperature * fabs(log(r1)));

      Point vn;
      vn = unit_normal * r;

      Point vni;
      vni = unit_normal * sp;

      Point vt;
      vt = vel - vni;
      double lvt = vt.Mod();
      Point tangent1;
      if (lvt < kEpsilon) {
        if (std::abs(unit_normal.x) < 0.9) {
          tangent1 = unit_normal.Cross(Point(1.0, 0.0, 0.0));
        } else {
          tangent1 = unit_normal.Cross(Point(0.0, 1.0, 0.0));
        }
        double tangent1_length = tangent1.Mod();
        if (tangent1_length < kEpsilon) {
          tangent1 = unit_normal.Cross(Point(0.0, 0.0, 1.0));
          tangent1_length = tangent1.Mod();
          if (tangent1_length < kEpsilon) {
            tangent1 = Point(0.0, 1.0, 0.0);
            tangent1_length = tangent1.Mod();
          }
        }
        tangent1 /= tangent1_length;
      } else {
        vt /= lvt;
        tangent1 = vt;
      }

      r = sqrt(2. * wall_temperature * fabs(log(r2)));
      double teta = 2. * Pi * r3;
      double vt1m = r * cos(teta);
      double vt2m = r * sin(teta);

      Point vt1 = tangent1;
      Point vt2;

      vt1 *= vt1m;

      vt2 = unit_normal.Cross(tangent1);
      double tangent2_length = vt2.Mod();
      if (tangent2_length < kEpsilon) {
        if (std::abs(unit_normal.z) < 0.9) {
          vt2 = unit_normal.Cross(Point(0.0, 0.0, 1.0));
        } else {
          vt2 = unit_normal.Cross(Point(0.0, 1.0, 0.0));
        }
        tangent2_length = vt2.Mod();
        if (tangent2_length < kEpsilon) {
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
                                    incoming_velocity - particle.velocity,
                                    0.5 * (incoming_velocity * incoming_velocity -
                                           particle.velocity * particle.velocity));

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

double InnerBoundary::CalcCellVolume(Point apex, Point size,
                                     utils::RandomStream&, Point* center) {
  const auto result = geometry_ ? geometry_->GasVolume(apex, size)
                                : CellVolume{size.Volume(), size.Volume(), 0,
                                             CellRegion::Gas, apex + size / 2};
  if (center) *center = result.gas_center.value_or(apex + size / 2);
  return result.gas;
}

void InnerBoundary::SetGeometry(Geometry& geometry) {
  geometry_ = &geometry;
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
  return geometry_;
}

bool InnerBoundary::Empty() const { return poligon_ptrs.empty(); }
}  // namespace mc3d
