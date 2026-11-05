#pragma once
#include "cell_cluster.h"
#include "config.h"
namespace mc3d {
class StatisticalDiagnostics {
 public:
  StatisticalDiagnostics(CellCluster& cluster, const SimulationConfig& config);
  void BeforeStep();
  void AfterStep();
  void Write();

 private:
  struct SurfaceSum {
    Point impulse{0, 0, 0};
    double heat = 0;
  };
  struct Window {
    double start, end, elapsed = 0;
    std::vector<WindowMoments> cells;
    std::vector<SurfaceSum> surface;
  };
  CellCluster& cluster_;
  bool statistics_, balance_;
  double fixed_temperature_ = 1, wall_initial_ = 0, max_wall_residual_ = 0;
  double weight_, start_ = 0, dt_ = 0, wall_before_ = 0;
  std::size_t steps_ = 0;
  ConservedQuantity before_, initial_, last_, incoming_, outgoing_;
  Point body_impulse_{0, 0, 0}, mirror_impulse_{0, 0, 0},
      momentum_residual_{0, 0, 0};
  double body_heat_ = 0, mirror_heat_ = 0, energy_residual_ = 0,
         max_energy_residual_ = 0, max_momentum_residual_ = 0;
  std::int64_t count_residual_ = 0;
  std::vector<Window> windows_;
  ConservedQuantity Inventory() const;
  double WallEnergy() const;
};
}  // namespace mc3d
