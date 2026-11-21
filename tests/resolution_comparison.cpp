#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

#include "cell_cluster.h"
#include "solver_setup.h"
namespace mc3d {
class CellClusterTestAccess {
 public:
  static auto& Cells(CellCluster& c) { return c.cells; }
};
class CellVolumeTestAccess {
 public:
  static auto& Particles(Cell& c) { return c.particles; }
};
}  // namespace mc3d
int main(int argc, char** argv) {
  if (argc != 2) return 2;
  std::ofstream out(argv[1]);
  out.precision(17);
  out << "seed;courant;time;xx;yy;zz;steps\n";
  std::array<std::array<double, 9>, 10> differences{};
  for (int trial = 0; trial < 10; ++trial)
    for (double courant : {.1, .05}) {
      mc3d::SimulationConfig config;
      mc3d::CellCluster cluster(1, 123 + trial);
      cluster.collision_courant = courant;
      cluster.collision_parameters = {
          mc3d::CollisionKind::VHS, 6.633521472521616e-26,
          3.8774087918641243e-10,   300,
          .7115290294715864,        300};
      cluster.SetApex({-.5, -.5, -.5});
      cluster.SetSize(1, 1, 1);
      cluster.SetEndTime(.1);
      cluster.Initialize(1, 1, 1, 200, .1, .3,
                         std::make_unique<mc3d::Geometry>(), 0, 0, 1);
      auto& cell = mc3d::CellClusterTestAccess::Cells(cluster).front();
      auto& particles = mc3d::CellVolumeTestAccess::Particles(cell);
      for (size_t i = 0; i < particles.size(); ++i)
        particles[i].velocity = {i % 2 ? std::sqrt(3.) : -std::sqrt(3.), 0, 0};
      for (int face = 0; face < 6; ++face)
        cluster.SetBoundaryCondition(mc3d::MakeBoundary(
            mc3d::MakeBoundaryDescriptor(static_cast<mc3d::BoundaryFace>(face),
                                         config),
            mc3d::BoundaryType::Periodic, config));
      int sample = 0;
      size_t steps = 0;
      for (double time : {.025, .05, .1}) {
        cluster.SetEndTime(time);
        bool finished = false;
        while (!finished) {
          finished = cluster.TimeStep().finished;
          if (++steps > 100000) return 3;
        }
        auto u = cell.ComputeMoments().velocity;
        std::array<double, 3> moments{};
        for (const auto& p : particles) {
          auto c = p.velocity - u;
          moments[0] += c.x * c.x;
          moments[1] += c.y * c.y;
          moments[2] += c.z * c.z;
        }
        out << 123 + trial << ';' << courant << ';' << time;
        for (int axis = 0; axis < 3; ++axis) {
          double m = moments[axis] / particles.size();
          out << ';' << m;
          differences[trial][sample * 3 + axis] += (courant == .1 ? 1 : -1) * m;
        }
        out << ';' << steps << '\n';
        ++sample;
      }
    }
  out.close();
  if (!out) return 4;
  for (int metric = 0; metric < 9; ++metric) {
    double mean = 0;
    for (auto& d : differences) mean += d[metric] / 10;
    double variance = 0;
    for (auto& d : differences)
      variance += (d[metric] - mean) * (d[metric] - mean) / 9;
    double tolerance = std::max(.06, 5 * std::sqrt(variance / 10));
    std::cout << "metric=" << metric << " difference=" << mean
              << " tolerance=" << tolerance << '\n';
    if (std::abs(mean) > tolerance) return 5;
  }
}
