#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "cell_cluster.h"
#include "config.h"
#include "free_boundary.h"
#include "geometry.h"
#include "giper_free_boundary.h"
#include "mirror_boundary.h"
#include "pereodic_boundary.h"
#include "point.h"

namespace {

struct BoundaryDescriptor {
  mc3d::Point position;
  mc3d::Point normal;
  mc3d::Point mixing;
};

BoundaryDescriptor MakeBoundaryDescriptor(mc3d::BoundaryFace face,
                                          const mc3d::SimulationConfig& cfg) {
  switch (face) {
    case mc3d::BoundaryFace::XNeg:
      return {mc3d::Point(-0.5 * cfg.Lx, 0.0, 0.0), mc3d::Point(-1.0, 0.0, 0.0),
              mc3d::Point(cfg.Ly, 0.0, 0.0)};
    case mc3d::BoundaryFace::XPos:
      return {mc3d::Point(0.5 * cfg.Lx, 0.0, 0.0), mc3d::Point(1.0, 0.0, 0.0),
              mc3d::Point(-cfg.Ly, 0.0, 0.0)};
    case mc3d::BoundaryFace::YNeg:
      return {mc3d::Point(0.0, -0.5 * cfg.Ly, 0.0), mc3d::Point(0.0, -1.0, 0.0),
              mc3d::Point(0.0, cfg.Lz, 0.0)};
    case mc3d::BoundaryFace::YPos:
      return {mc3d::Point(0.0, 0.5 * cfg.Ly, 0.0), mc3d::Point(0.0, 1.0, 0.0),
              mc3d::Point(0.0, -cfg.Lz, 0.0)};
    case mc3d::BoundaryFace::ZNeg:
      return {mc3d::Point(0.0, 0.0, -0.5 * cfg.Lz), mc3d::Point(0.0, 0.0, -1.0),
              mc3d::Point(0.0, 0.0, cfg.Lz)};
    case mc3d::BoundaryFace::ZPos:
      return {mc3d::Point(0.0, 0.0, 0.5 * cfg.Lz), mc3d::Point(0.0, 0.0, 1.0),
              mc3d::Point(0.0, 0.0, -cfg.Lz)};
  }
  return {};
}

std::unique_ptr<mc3d::Boundary> MakeBoundary(
    const BoundaryDescriptor& descriptor, mc3d::BoundaryType type,
    const mc3d::SimulationConfig& cfg) {
  switch (type) {
    case mc3d::BoundaryType::None:
      return nullptr;
    case mc3d::BoundaryType::Mirror:
      return std::make_unique<mc3d::MirrorBoundary>(descriptor.position,
                                                    descriptor.normal);
    case mc3d::BoundaryType::Periodic:
      return std::make_unique<mc3d::PeriodicBoundary>(
          descriptor.position, descriptor.normal, descriptor.mixing);
    case mc3d::BoundaryType::Free:
      return std::make_unique<mc3d::FreeBoundary>(
          descriptor.position, descriptor.normal, cfg.particles_per_cell, cfg.S,
          cfg.temperature, cfg.alpha);
    case mc3d::BoundaryType::HyperFree:
      return std::make_unique<mc3d::HyperFreeBoundary>(
          descriptor.position, descriptor.normal, cfg.particles_per_cell, cfg.S,
          cfg.temperature, cfg.alpha);
  }
  return nullptr;
}

void ApplyGeometryPostProcessing(const mc3d::SimulationConfig& cfg,
                                 mc3d::Geometry& geometry) {
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
    const mc3d::Point delta(cfg.geometry_move_x.value_or(0.0),
                            cfg.geometry_move_y.value_or(0.0),
                            cfg.geometry_move_z.value_or(0.0));
    geometry.Move(delta);
  }
}

std::unique_ptr<mc3d::Geometry> BuildGeometry(
    const mc3d::SimulationConfig& cfg) {
  if (!cfg.geometry_file.empty()) {
    auto body = std::make_unique<mc3d::Geometry>(cfg.geometry_file.c_str());
    ApplyGeometryPostProcessing(cfg, *body);
    return body;
  }

  auto body = std::make_unique<mc3d::Geometry>();

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
  }

  ApplyGeometryPostProcessing(cfg, *body);
  return body;
}

}  // namespace

int main(int argc, char* argv[]) {
  const auto config = mc3d::LoadSimulationConfig(argc, argv);

  auto geometry = BuildGeometry(config);

  const mc3d::Point apex(config.apex_x.value_or(-0.5 * config.Lx),
                         config.apex_y.value_or(-0.5 * config.Ly),
                         config.apex_z.value_or(-0.5 * config.Lz));

  mc3d::CellCluster cluster;
  cluster.SetApex(apex);
  cluster.SetSize(config.Lx, config.Ly, config.Lz);

  const double total_particles =
      static_cast<double>(config.particles_per_cell) *
      static_cast<double>(config.cells_x) *
      static_cast<double>(config.cells_y) * static_cast<double>(config.cells_z);

  cluster.Initialize(config.cells_x, config.cells_y, config.cells_z,
                     total_particles, config.Kn, config.Cu, std::move(geometry),
                     config.S, config.alpha, config.temperature);
  cluster.SetBinaryOutput(config.snapshots_binary);
  cluster.SetSnapshotInterval(config.snapshot_interval);

  if (config.write_default_snapshot_before_compute) {
    cluster.WriteFile();
  }

  std::vector<std::unique_ptr<mc3d::Boundary>> boundaries;
  boundaries.reserve(config.boundary_types.size());

  for (std::size_t idx = 0; idx < config.boundary_types.size(); ++idx) {
    const auto face = static_cast<mc3d::BoundaryFace>(idx);
    const auto descriptor = MakeBoundaryDescriptor(face, config);
    auto boundary =
        MakeBoundary(descriptor, config.boundary_types[idx], config);
    if (boundary) {
      boundaries.emplace_back(std::move(boundary));
    }
  }

  if (!boundaries.empty()) {
    cluster.SetBoundaryCondition(std::move(boundaries));
  }

  if (config.write_speed_before) {
    cluster.WriteSpeedFile();
  }

  cluster.SetEndTime(config.end_time);

  if (!config.precompute_snapshot.empty()) {
    cluster.WriteFile(config.precompute_snapshot);
  }

  cluster.Compute();

  std::cout << "Computation is over." << std::endl;

  if (!config.final_cell_snapshot.empty()) {
    cluster.WriteCellFile(config.final_cell_snapshot);
  }

  if (config.write_default_snapshot_after_compute) {
    cluster.WriteFile();
  }

  if (!config.final_snapshot.empty()) {
    cluster.WriteFile(config.final_snapshot);
  }

  if (config.write_times) {
    cluster.WriteTimes();
  }

  if (config.write_speed_after) {
    cluster.WriteSpeedFile();
  }

  return 0;
}