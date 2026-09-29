#include "solver_setup.h"

#include <utility>
#include <cmath>
#include <stdexcept>

#include "free_boundary.h"
#include "giper_free_boundary.h"
#include "mirror_boundary.h"
#include "pereodic_boundary.h"

namespace mc3d {

double ReferenceParticleDensity(const SimulationConfig& cfg) {
  const double volume = cfg.Lx * cfg.Ly * cfg.Lz;
  if (!(volume > 0) || !std::isfinite(volume))
    throw std::invalid_argument("Domain volume must be finite and positive");
  return static_cast<double>(cfg.particles_per_cell) * cfg.cells_x *
         cfg.cells_y * cfg.cells_z / volume;
}

BoundaryDescriptor MakeBoundaryDescriptor(BoundaryFace face,
                                          const SimulationConfig& cfg) {
  const Point apex(cfg.apex_x.value_or(-0.5 * cfg.Lx),
                   cfg.apex_y.value_or(-0.5 * cfg.Ly),
                   cfg.apex_z.value_or(-0.5 * cfg.Lz));
  switch (face) {
    case BoundaryFace::XNeg:
      return {Point(apex.x, 0, 0), Point(-1, 0, 0),
              Point(cfg.Lx, 0, 0)};
    case BoundaryFace::XPos:
      return {Point(apex.x + cfg.Lx, 0, 0), Point(1, 0, 0),
              Point(-cfg.Lx, 0, 0)};
    case BoundaryFace::YNeg:
      return {Point(0, apex.y, 0), Point(0, -1, 0),
              Point(0, cfg.Ly, 0)};
    case BoundaryFace::YPos:
      return {Point(0, apex.y + cfg.Ly, 0), Point(0, 1, 0),
              Point(0, -cfg.Ly, 0)};
    case BoundaryFace::ZNeg:
      return {Point(0, 0, apex.z), Point(0, 0, -1),
              Point(0, 0, cfg.Lz)};
    case BoundaryFace::ZPos:
      return {Point(0, 0, apex.z + cfg.Lz), Point(0, 0, 1),
              Point(0, 0, -cfg.Lz)};
  }
  return {};
}

std::unique_ptr<Boundary> MakeBoundary(const BoundaryDescriptor& descriptor,
                                       BoundaryType type,
                                       const SimulationConfig& cfg) {
  switch (type) {
    case BoundaryType::None:
      return nullptr;
    case BoundaryType::Mirror:
      return std::make_unique<MirrorBoundary>(descriptor.position,
                                              descriptor.normal);
    case BoundaryType::Periodic:
      return std::make_unique<PeriodicBoundary>(
          descriptor.position, descriptor.normal, descriptor.mixing);
    case BoundaryType::Free:
      return std::make_unique<FreeBoundary>(
          descriptor.position, descriptor.normal, cfg.particles_per_cell, cfg.S,
          cfg.temperature, cfg.alpha);
    case BoundaryType::HyperFree:
      return std::make_unique<HyperFreeBoundary>(
          descriptor.position, descriptor.normal, cfg.particles_per_cell, cfg.S,
          cfg.temperature, cfg.alpha);
  }
  return nullptr;
}

void ApplyGeometryPostProcessing(const SimulationConfig& cfg,
                                 Geometry& geometry) {
  if (cfg.geometry_fix_polygons) {
    geometry.FixPolygons();
  }
  if (cfg.geometry_reverse_normals) {
    geometry.ReverseNormals();
  }
  if (cfg.geometry_scale) {
    geometry.Scale(*cfg.geometry_scale);
  }
  if (cfg.geometry_fragment_length) {
    geometry.Fragment(*cfg.geometry_fragment_length);
  }
  if (cfg.geometry_move_x || cfg.geometry_move_y || cfg.geometry_move_z) {
    const Point delta(cfg.geometry_move_x.value_or(0.0),
                      cfg.geometry_move_y.value_or(0.0),
                      cfg.geometry_move_z.value_or(0.0));
    geometry.Move(delta);
  }
}

std::unique_ptr<Geometry> BuildGeometry(const SimulationConfig& cfg) {
  if (!cfg.geometry_file.empty()) {
    auto body = std::make_unique<Geometry>(cfg.geometry_file.c_str());
    ApplyGeometryPostProcessing(cfg, *body);
    return body;
  }

  auto body = std::make_unique<Geometry>();

  if (cfg.geometry_type == "wedge") {
    body->CreateWedge(cfg.geometry_wedge_x, cfg.geometry_wedge_width,
                      cfg.geometry_wedge_length, cfg.geometry_wedge_alpha);
  } else if (cfg.geometry_type == "pyramid") {
    body->CreatePyramid(cfg.geometry_pyramid_x, cfg.geometry_pyramid_width,
                        cfg.geometry_pyramid_length,
                        cfg.geometry_pyramid_height);
  } else if (cfg.geometry_type == "cube") {
    body->CreateCube(cfg.geometry_cube_x, cfg.geometry_cube_width,
                     cfg.geometry_cube_length, cfg.geometry_cube_height);
  } else if (cfg.geometry_type == "cylinder") {
    body->CreateCylinder(cfg.geometry_cylinder_x, cfg.geometry_cylinder_radius,
                         cfg.geometry_cylinder_length,
                         cfg.geometry_cylinder_segments);
  }

  ApplyGeometryPostProcessing(cfg, *body);
  return body;
}

}  // namespace mc3d
