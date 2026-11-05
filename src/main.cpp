#include "simulation_runner.h"
#include "run_manifest.h"
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

#include "cell_cluster.h"
#include "config.h"
#include "point.h"
#include "solver_setup.h"

int Run(int argc, char* argv[]) {
  auto config = mc3d::LoadSimulationConfig(argc, argv);
  if (!config.seed) config.seed = std::random_device{}();
  config.thread_pool_size = std::max(1u, config.thread_pool_size);
  mc3d::RunManifest manifest(config, argc, argv);
  std::signal(SIGINT, [](int signal) { mc3d::stop_requested = signal; });
  std::signal(SIGTERM, [](int signal) { mc3d::stop_requested = signal; });
  try {

  auto geometry = mc3d::BuildGeometry(config);
  manifest.GeometryHash(*geometry);

  const mc3d::Point apex(config.apex_x.value_or(-0.5 * config.Lx),
                         config.apex_y.value_or(-0.5 * config.Ly),
                         config.apex_z.value_or(-0.5 * config.Lz));

  const auto seed = config.seed ? *config.seed : std::random_device{}();
  std::cout << "Random generator: " << utils::RandomStream::kAlgorithm << std::endl;
  std::cout << "Random seed: " << seed << std::endl;
  mc3d::CellCluster cluster(config.thread_pool_size, seed);
  mc3d::SimulationRunner runner(cluster);
  cluster.collision_parameters = config.collisions;
  cluster.SetApex(apex);
  cluster.SetSize(config.Lx, config.Ly, config.Lz);

  const double reference_density = mc3d::ReferenceParticleDensity(config);

  cluster.Initialize(config.cells_x, config.cells_y, config.cells_z,
                     reference_density, config.Kn, config.Cu, std::move(geometry),
                     config.S, config.alpha, config.temperature, config.wall_temperature,
                     config.wall_heat_capacity);
  if (config.force_reference_area)
    runner.SetForceReference(*config.force_reference_area,
                             config.S * std::sqrt(2 * config.temperature), config.alpha);
  runner.Writer().SetBinaryOutput(config.snapshots_binary);
  runner.SetSnapshotInterval(config.snapshot_interval);

  if (config.write_default_snapshot_before_compute) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteFile());
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
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteSpeedFile());
  }

  cluster.SetEndTime(config.end_time);

  if (!config.precompute_snapshot.empty()) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteFile(config.precompute_snapshot));
  }

  runner.progress = [&manifest](double t, std::size_t steps) { manifest.Progress(t, steps); };
  runner.runtime_progress = [&manifest](double t, std::size_t step, const std::string& mode, std::size_t workers) { manifest.RuntimePolicy(t, step, mode, workers); };
  runner.SetRuntimeControl(config.runtime_control);
  runner.Compute();
  manifest.Progress(runner.Time(), runner.Steps());

  std::cout << "Computation is over." << std::endl;

  if (!config.final_cell_snapshot.empty()) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteCellFile(config.final_cell_snapshot));
  }

  if (config.write_default_snapshot_after_compute) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteFile());
  }

  if (!config.final_snapshot.empty()) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteFile(config.final_snapshot));
  }

  if (config.write_times) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteTimes());
  }

  if (config.write_speed_after) {
    mc3d::SimulationRunner::RequireWrite(runner.Writer().WriteSpeedFile());
  }

  manifest.Finish(mc3d::stop_requested ? "interrupted" : "completed");
  return mc3d::stop_requested ? 128 + mc3d::stop_requested : 0;
  } catch(const std::exception& error) { manifest.Finish("failed", error.what()); throw; }
}

int main(int argc, char* argv[]) {
  try { return Run(argc, argv); }
  catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
