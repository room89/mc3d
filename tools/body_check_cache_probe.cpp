// Build against current mc3d_core; reports cache metadata without instrumenting solver logs.
#include "cell_cluster.h"
#include "config.h"
#include "solver_setup.h"
#include <iostream>
namespace mc3d {
class CellClusterTestAccess {
 public:
  static void Probe(CellCluster& c, const SimulationConfig& cfg) {
    c.body_ = BuildGeometry(cfg);
    if (c.body_) c.body_->PrepareQueries();
    const Point origin(cfg.apex_x.value_or(-.5*cfg.Lx), cfg.apex_y.value_or(-.5*cfg.Ly), cfg.apex_z.value_or(-.5*cfg.Lz));
    const Point size(cfg.Lx/cfg.cells_x, cfg.Ly/cfg.cells_y, cfg.Lz/cfg.cells_z);
    for (unsigned i=0;i<cfg.cells_x;++i) for (unsigned j=0;j<cfg.cells_y;++j) for (unsigned k=0;k<cfg.cells_z;++k) {
      c.cells.emplace_back(123,c.cells.size());
      c.cells.back().SetSize(size);
      c.cells.back().SetApex(origin+Point(i*size.x,j*size.y,k*size.z));
    }
    c.CountBodyParticles();
    std::cout << "cells=" << c.cells.size() << " candidates=" << c.body_check_cells_.size()
              << " cache_capacity_bytes=" << c.body_check_cells_.capacity()*sizeof(Cell*) << '\n';
  }
};
}
int main(int argc, char** argv) {
  const auto cfg = mc3d::LoadSimulationConfig(argc, argv);
  mc3d::CellCluster cluster(1,123);
  mc3d::CellClusterTestAccess::Probe(cluster,cfg);
}
