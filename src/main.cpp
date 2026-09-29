#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "cell_cluster.h"
#include "config.h"
#include "point.h"
#include "solver_setup.h"

int main(int argc, char* argv[]) {
  const auto config = mc3d::LoadSimulationConfig(argc, argv);

  auto geometry = mc3d::BuildGeometry(config);

  const mc3d::Point apex(config.apex_x.value_or(-0.5 * config.Lx),
                         config.apex_y.value_or(-0.5 * config.Ly),
                         config.apex_z.value_or(-0.5 * config.Lz));

  mc3d::CellCluster cluster(config.thread_pool_size);
  cluster.SetApex(apex);
  cluster.SetSize(config.Lx, config.Ly, config.Lz);

  const double reference_density = mc3d::ReferenceParticleDensity(config);

  cluster.Initialize(config.cells_x, config.cells_y, config.cells_z,
                     reference_density, config.Kn, config.Cu, std::move(geometry),
                     config.S, config.alpha, config.temperature, config.wall_temperature);
  cluster.SetBinaryOutput(config.snapshots_binary);
  cluster.SetSnapshotInterval(config.snapshot_interval);

  if (config.write_default_snapshot_before_compute) {
    cluster.WriteFile();
  }

  std::vector<std::unique_ptr<mc3d::Boundary>> boundaries;
  boundaries.reserve(config.boundary_types.size());

  for (std::size_t idx = 0; idx < config.boundary_types.size(); ++idx) {
    const auto face = static_cast<mc3d::BoundaryFace>(idx);
    const auto descriptor = mc3d::MakeBoundaryDescriptor(face, config);
    auto boundary =
        mc3d::MakeBoundary(descriptor, config.boundary_types[idx], config);
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
