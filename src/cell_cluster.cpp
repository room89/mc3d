#include "cell_cluster.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <exception>
#include <future>
#include <iterator>
#include <sstream>
#include <string>
#include <utils/logger.hpp>
#include <utils/stopwatch.hpp>

using namespace std;
using namespace mc3d;

namespace mc3d {
CellCluster::CellCluster(std::size_t thread_pool_size)
    : thread_pool_(std::max<std::size_t>(1, thread_pool_size),
                   "move_and_collisions") {
  Kn = 0;
  np = 0;
  ncx = 0;
  ncy = 0;
  ncz = 0;
  N = 0;
  Lx = 0;
  Ly = 0;
  Lz = 0;
  t = 0;
  dt = 10000;
  step = 0;
}

CellCluster::~CellCluster() {
  thread_pool_.Stop();
  cells.clear();
}

bool CellCluster::Initialize(unsigned int ncx, unsigned int ncy,
                             unsigned int ncz, double density, double Kn,
                             double Cu, std::unique_ptr<Geometry>&& body,
                             double S, double alpha, double T, double wall_temperature) {
  LOG_DEBUG() << bool(body);
  body_ = std::move(body);

  srand(static_cast<unsigned int>(time(nullptr)));

  this->ncx = ncx;
  this->ncy = ncy;
  this->ncz = ncz;

  this->density = density;

  this->Kn = Kn;
  this->Cu = Cu;

  unsigned int N = 0;

  double dx = Lx / double(ncx);
  double dy = Ly / double(ncy);
  double dz = Lz / double(ncz);
  cell_dx_ = dx;
  cell_dy_ = dy;
  cell_dz_ = dz;
  cell_lookup_.clear();

  this->t = 0;
  this->dt = 1000000000;

  const size_t total_cells = static_cast<size_t>(ncx) *
                             static_cast<size_t>(ncy) *
                             static_cast<size_t>(ncz);
  size_t initialized_cells = 0;
  size_t next_progress_percent = 1;
  if (total_cells > 0) {
    LOG_INFO() << "Initializing " << total_cells << " cells";
  }

  auto cell_iter = cells.begin();
  for (size_t i = 0; i < ncx; i++) {
    for (size_t j = 0; j < ncy; j++) {
      for (size_t k = 0; k < ncz; k++) {
        Cell temp_cell;

        size_t N_ = 0;
        Point a;
        a.x = apex.x + i * dx;
        a.y = apex.y + j * dy;
        a.z = apex.z + k * dz;

        N_ = static_cast<unsigned int>(density * dx * dy * dz);
        temp_cell.SetParameters(S, alpha, T);
        temp_cell.SetWallTemperature(wall_temperature);

        temp_cell.SetApex(a);
        temp_cell.SetSize(dx, dy, dz);

        temp_cell.Initialize(N_, body_);
        if (density > 0)
          temp_cell.SetReferenceParticleCount(density * dx * dy * dz);

        temp_cell.SetCharacteristicLength(Lx);
        temp_cell.SetKn(Kn);

        N += N_;
        double dtt = temp_cell.GetDt();
        dt = min(dt, dtt);

        cells.emplace_back(std::move(temp_cell));
        ++initialized_cells;
        if (total_cells > 0) {
          size_t completed_percent = (initialized_cells * 100) / total_cells;
          while (next_progress_percent <= 100 &&
                 completed_percent >= next_progress_percent) {
            LOG_INFO() << "Cell initialization progress: "
                       << next_progress_percent << "% (" << initialized_cells
                       << "/" << total_cells << ")";
            ++next_progress_percent;
          }
        }
      }
    }
  }

  //  cells_test();

  this->N = N;

  cell_lookup_.clear();
  cell_lookup_.reserve(cells.size());
  for (auto& cell : cells) {
    cell_lookup_.push_back(&cell);
  }

  // density is the reservoir number of model particles per unit volume.
  // It must not become the total population (nor the population outside a body).

  cell_iter = cells.begin();
  for (auto& cell : cells) {
    cell.CleanInnerParticles(*body_);
  }

  this->SyncDt();
  cout << "Cluster initialize is done" << endl;

  return true;
}

void CellCluster::SetApex(Point apex) { this->apex = apex; }

void CellCluster::SetSize(double Lx, double Ly, double Lz) {
  this->Lx = Lx;
  this->Ly = Ly;
  this->Lz = Lz;
}

StepResult CellCluster::TimeStep() {
  SyncDt();

  LOG_INFO() << "dt = " << dt << " t = " << t;

  //  auto cell_iter = cells.begin();

  {
    utils::Stopwatch sw("TimeStep::CalculateCells");
    std::vector<std::future<void>> futures;
    futures.reserve(cells.size());
    for (auto& cell : cells) {
      futures.emplace_back(thread_pool_.Execute(
          [cell_ptr = &cell]() { cell_ptr->Calculate(); }));
    }
    // A failing cell must not leave other tasks using cells after the caller
    // catches the error or destroys the cluster.
    std::exception_ptr failure;
    for (auto& fut : futures) {
      try {
        fut.get();
      } catch (...) {
        if (!failure) failure = std::current_exception();
      }
    }
    if (failure) std::rethrow_exception(failure);
  }

  {
    utils::Stopwatch sw("TimeStep::GatherParticles");
    for (auto& cell : cells) {
      cell.Sort();
      cell.SortNeighbors();
      auto& cell_buffer = cell.GetBuffer();
      partile_buffer.insert(partile_buffer.end(),
                            std::make_move_iterator(cell_buffer.begin()),
                            std::make_move_iterator(cell_buffer.end()));
      cell_buffer.clear();
    }
    LOG_INFO() << "particle_buffer_size:" << partile_buffer.size();
  }

  {
    utils::Stopwatch sw("TimeStep::AddParticlesBeforeBoundary");
    DistributeParticles(partile_buffer);
  }

  {
    utils::Stopwatch sw("TimeStep::BoundaryConditions");
    BoundaryCondition();
  }

  {
    utils::Stopwatch sw("TimeStep::AddParticlesAfterBoundary");
    DistributeParticles(partile_buffer);
  }
  t += dt;

  partile_buffer.clear();

  std::size_t inner_particles = 0;
  if (body_ && body_->PolygonCount() > 0) {
    const auto [lower, upper] = body_->Bounds();
    for (const auto& cell : cells) {
      const auto origin = cell.GetApex();
      const auto extent = cell.GetSize();
      if (origin.x <= upper.x && origin.x + extent.x >= lower.x &&
          origin.y <= upper.y && origin.y + extent.y >= lower.y &&
          origin.z <= upper.z && origin.z + extent.z >= lower.z) {
        inner_particles += cell.CountInnerParticles(*body_);
      }
    }
  }
  return {t, dt, t >= t_end, inner_particles};
}

void CellCluster::BoundaryCondition() {
  for (auto& boundary : boundary_cond_outer) {
    boundary->BoundaryCondition(partile_buffer, dt);
  }
}

bool CellCluster::FindCellIndex(const Point& position,
                                size_t& cell_index) const {
  if (cell_lookup_.empty() || cell_dx_ <= 0.0 || cell_dy_ <= 0.0 ||
      cell_dz_ <= 0.0 || ncx == 0 || ncy == 0 || ncz == 0) {
    return false;
  }

  if (position.x < apex.x || position.x >= apex.x + Lx ||
      position.y < apex.y || position.y >= apex.y + Ly ||
      position.z < apex.z || position.z >= apex.z + Lz) {
    return false;
  }

  const double rel_x = (position.x - apex.x) / cell_dx_;
  const double rel_y = (position.y - apex.y) / cell_dy_;
  const double rel_z = (position.z - apex.z) / cell_dz_;

  unsigned int ix = static_cast<unsigned int>(rel_x);
  unsigned int iy = static_cast<unsigned int>(rel_y);
  unsigned int iz = static_cast<unsigned int>(rel_z);

  if (ix >= ncx) ix = ncx - 1;
  if (iy >= ncy) iy = ncy - 1;
  if (iz >= ncz) iz = ncz - 1;

  cell_index = (static_cast<size_t>(ix) * ncy + iy) * ncz + iz;
  return cell_index < cell_lookup_.size();
}

void CellCluster::DistributeParticles(std::vector<Particle>& buffer) {
  if (buffer.empty()) {
    return;
  }

  const size_t thread_count =
      std::max<size_t>(1, thread_pool_.GetThreadCount());
  const size_t chunk_size =
      std::max<size_t>(1024, (buffer.size() + thread_count - 1) / thread_count);
  const size_t chunk_count = (buffer.size() + chunk_size - 1) / chunk_size;

  std::vector<std::future<void>> futures;
  futures.reserve(chunk_count);
  std::vector<std::vector<Particle>> leftovers(chunk_count);

  size_t chunk_idx = 0;
  for (size_t start = 0; start < buffer.size();
       start += chunk_size, ++chunk_idx) {
    const size_t end = std::min(buffer.size(), start + chunk_size);
    auto& chunk_leftovers = leftovers[chunk_idx];
    chunk_leftovers.reserve(end - start);
    futures.emplace_back(
        thread_pool_.Execute([this, &buffer, start, end, &chunk_leftovers]() {
          for (size_t idx = start; idx < end; ++idx) {
            Particle particle = std::move(buffer[idx]);
            size_t cell_index = 0;
            bool accepted = false;
            if (FindCellIndex(particle.position, cell_index)) {
              Cell* cell_ptr = cell_lookup_[cell_index];
              if (cell_ptr != nullptr) {
                accepted = cell_ptr->TryAcceptParticle(particle);
              }
            }
            if (!accepted) {
              chunk_leftovers.push_back(std::move(particle));
            }
          }
        }));
  }

  for (auto& future : futures) {
    future.get();
  }

  buffer.clear();
  size_t total_remaining = 0;
  for (auto& chunk_leftovers : leftovers) {
    total_remaining += chunk_leftovers.size();
  }
  buffer.reserve(total_remaining);
  for (auto& chunk_leftovers : leftovers) {
    buffer.insert(buffer.end(),
                  std::make_move_iterator(chunk_leftovers.begin()),
                  std::make_move_iterator(chunk_leftovers.end()));
  }
}

bool CellCluster::SendData() { return true; }

bool CellCluster::RecvData() { return true; }

void CellCluster::SetBoundaryCondition(
    std::vector<std::unique_ptr<Boundary>>&& boundaries) {
  for (auto& boundary : boundaries) {
    SetBoundaryCondition(std::move(boundary));
  }
}

void CellCluster::SetBoundaryCondition(std::unique_ptr<Boundary> boundary) {
  if (!boundary) {
    return;
  }

  if (auto* free_boundary = dynamic_cast<FreeBoundary*>(boundary.get())) {
    free_boundary->AddCell(cells);
  }

  boundary_cond_outer.push_back(std::move(boundary));
}

void CellCluster::SetEndTime(double t_end) {
  this->data_t = 0;
  this->data_dt = t_end / 100 + 0.000001;
  this->t_end = t_end;
}
void CellCluster::CalculateDt() {
  dt = 1000000.;

  cell_iter = cells.begin();
  for (auto& cell : cells) {
    double dtt = cell.CalculateDt();
    if (dtt < dt) dt = dtt;
  }
  dt *= Cu;
  if (t + dt > t_end) dt = t_end - t;
}

void CellCluster::SetDtInCells(double dt) {
  for (auto& cell : cells) {
    cell.SetDt(dt);
  }
}

void CellCluster::SyncDt() {
  CalculateDt();
  SetDtInCells(dt);
}

void CellCluster::SyncData() {}

void CellCluster::FragmentCells() {}

void CellCluster::TestCells() {
  cell_iter = cells.begin();
  while (cell_iter != cells.end()) {
    cell_iter++;
  }
}

void CellCluster::SetDataSaveDtime(double data_dt) { this->data_dt = data_dt; }

void CellCluster::CleanInnerParticles() {
  for (auto& cell : cells) {
    cell.CleanInnerParticles(*body_);
  }
}
}  // namespace mc3d
