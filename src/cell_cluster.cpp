#include "cell_cluster.h"

#include <algorithm>
#include <atomic>
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
CellCluster::CellCluster(std::size_t thread_pool_size, uint32_t seed)
    : random_seed_(seed), thread_pool_(std::max<std::size_t>(1, thread_pool_size),
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
                             double S, double alpha, double T, double wall_temperature,
                             std::optional<double> wall_heat_capacity) {
  LOG_DEBUG() << bool(body);
  InvalidateBodyCheck();
  body_ = std::move(body);
  if (body_) body_->PrepareQueries();
  if (wall_heat_capacity) {
    if (!body_ || !(density > 0) || !std::isfinite(density))
      throw std::invalid_argument("Surface heating requires a body and positive reference particle density");
    body_->ConfigureSurfaceHeating(wall_temperature, *wall_heat_capacity, 1.0 / density);
  }

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
  body_impulse_ = {0, 0, 0};
  if (body_) body_->TakeStepImpulse();
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
        cells.emplace_back(random_seed_, initialized_cells);
        Cell& cell = cells.back();

        size_t N_ = 0;
        Point a;
        a.x = apex.x + i * dx;
        a.y = apex.y + j * dy;
        a.z = apex.z + k * dz;

        N_ = static_cast<unsigned int>(density * dx * dy * dz);
        cell.SetParameters(S, alpha, T);
        cell.SetWallTemperature(wall_temperature);

        cell.SetApex(a);
        cell.SetSize(dx, dy, dz);

        cell.Initialize(N_, body_);
        cell.defer_surface_exchange = true;
        if (density > 0)
          cell.SetReferenceParticleCount(density * dx * dy * dz);

        cell.SetCharacteristicLength(Lx);
        cell.SetKn(Kn);
        cell.collision_parameters = collision_parameters;

        N += N_;
        double dtt = cell.GetDt();
        dt = min(dt, dtt);

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
  this->N = 0;
  for (auto& cell : cells) {
    if (body_) cell.CleanInnerParticles(*body_);
    this->N += cell.GetParticleCount();
  }

  this->SyncDt();
  cout << "Cluster initialize is done" << endl;

  return true;
}

void CellCluster::SetApex(Point apex) {
  InvalidateBodyCheck();
  this->apex = apex;
}

void CellCluster::SetSize(double Lx, double Ly, double Lz) {
  InvalidateBodyCheck();
  this->Lx = Lx;
  this->Ly = Ly;
  this->Lz = Lz;
}

StepResult CellCluster::TimeStep() {
  if (body_) body_->ResetStepSurfaceExchange();
  step_incoming = {};
  step_outgoing = {};
  step_body_impulse = {0, 0, 0};
  step_mirror_impulse = {0, 0, 0};
  step_body_heat = step_mirror_heat = 0;
  for (auto& cell : cells) {
    cell.incoming_quantity = {};
    cell.surface_exchanges.clear();
  }
  for (auto& boundary : boundary_cond_outer) {
    boundary->step_impulse = {0, 0, 0};
    boundary->step_heat = 0;
  }
  {
    utils::Stopwatch sw("TimeStep::SyncDt");
    SyncDt();
  }

  if (t < t_end && (!(dt > 0) || !std::isfinite(dt) || t + dt == t))
    throw std::runtime_error("Time step cannot advance simulation time");
  if (before_step) before_step();
  LOG_INFO() << "dt = " << dt << " t = " << t;

  //  auto cell_iter = cells.begin();

  {
    utils::Stopwatch sw("TimeStep::CalculateCells");
    const auto workers = std::max<std::size_t>(1,
        std::min(thread_pool_.GetThreadCount(), cells.size()));
    // Keep several chunks per worker even on small grids; cap large-grid chunks.
    const auto chunk_size = std::clamp(
        (cells.size() + 4 * workers - 1) / (4 * workers),
        std::size_t{1}, std::size_t{32});
    std::atomic<std::size_t> next{0};
    auto calculate = [&] {
      std::exception_ptr failure;
      for (;;) {
        const auto begin = next.fetch_add(chunk_size, std::memory_order_relaxed);
        if (begin >= cells.size()) break;
        const auto end = std::min(begin + chunk_size, cells.size());
        for (auto i = begin; i < end; ++i) {
          // SyncDt prepared moments; only the common dt changed since then.
          try { cells[i].CalculateWithPreparedMoments(); }
          catch (...) { if (!failure) failure = std::current_exception(); }
        }
      }
      if (failure) std::rethrow_exception(failure);
    };
    if (workers <= 1) {
      calculate();
    } else {
      std::vector<std::future<void>> futures;
      futures.reserve(workers);
      std::exception_ptr failure;
      try {
        for (std::size_t i = 0; i < workers; ++i)
          futures.emplace_back(thread_pool_.Execute(calculate));
      } catch (...) { failure = std::current_exception(); }
      // Drain every submitted task, including when submission itself fails.
      for (auto& future : futures) {
        try { future.get(); }
        catch (...) { if (!failure) failure = std::current_exception(); }
      }
      if (failure) std::rethrow_exception(failure);
    }
  }

  OutgoingSources outgoing;
  {
    utils::Stopwatch sw("TimeStep::GatherParticles");
    outgoing = GatherParticles();
    LOG_INFO() << "particle_buffer_size:" << (outgoing.empty() ? 0 : outgoing.back().end);
  }

  {
    utils::Stopwatch sw("TimeStep::AddParticlesBeforeBoundary");
    DistributeSources(outgoing, partile_buffer);
    outgoing.clear();
  }

  {
    utils::Stopwatch sw("TimeStep::BoundaryConditions");
    BoundaryCondition();
  }

  {
    utils::Stopwatch sw("TimeStep::AddParticlesAfterBoundary");
    DistributeParticles(partile_buffer);
  }
  for (const auto& cell : cells) cell.ValidateGasPopulation();
  // Commit surface exchanges in cell/event order, independently of worker
  // scheduling.
  for (auto& cell : cells) {
    for (const auto& e : cell.surface_exchanges)
      e.body->AccumulateForce(*e.polygon, e.impulse, e.energy);
    cell.surface_exchanges.clear();
  }
  if (body_) body_->UpdateSurfaceTemperatures();
  t += dt;

  std::array<bool, 6> may_escape{true, true, true, true, true, true};
  for (const auto& boundary : boundary_cond_outer) {
    if (!dynamic_cast<PeriodicBoundary*>(boundary.get()) &&
        !dynamic_cast<MirrorBoundary*>(boundary.get()))
      continue;
    const auto normal = boundary->GetNormal();
    const int axis = normal.x != 0 ? 0 : normal.y != 0 ? 1 : 2;
    const double component = axis == 0   ? normal.x
                             : axis == 1 ? normal.y
                                         : normal.z;
    may_escape[2 * axis + (component > 0 ? 1 : 0)] = false;
  }
  for (const auto& particle : partile_buffer) {
    if (!std::isfinite(particle.position.x) ||
        !std::isfinite(particle.position.y) ||
        !std::isfinite(particle.position.z))
      throw std::runtime_error("Nonfinite outgoing particle");
    size_t outgoing_index = 0;
    if (FindCellIndex(particle.position, outgoing_index))
      throw std::runtime_error("Unrouted particle inside domain");
    const auto p = particle.position;
    const bool expected = (p.x < apex.x && may_escape[0]) ||
                          (p.x >= apex.x + Lx && may_escape[1]) ||
                          (p.y < apex.y && may_escape[2]) ||
                          (p.y >= apex.y + Ly && may_escape[3]) ||
                          (p.z < apex.z && may_escape[4]) ||
                          (p.z >= apex.z + Lz && may_escape[5]);
    if (!expected)
      throw std::runtime_error("Particle escaped a closed boundary");
    step_outgoing.Add(particle.velocity);
  }
  for (const auto& cell : cells) step_incoming += cell.incoming_quantity;
  for (const auto& boundary : boundary_cond_outer) {
    step_mirror_impulse += boundary->step_impulse;
    step_mirror_heat += boundary->step_heat;
  }
  partile_buffer.clear();

  std::size_t inner_particles = 0;
  if (body_ && body_->PolygonCount() > 0) {
    utils::Stopwatch sw("TimeStep::CheckBodyParticles");
    inner_particles = CountBodyParticles();
  } else {
    InvalidateBodyCheck();
  }
  Point force(0, 0, 0), mean_force(0, 0, 0);
  if (body_ && body_->PolygonCount() > 0) {
    step_body_impulse = body_->TakeStepImpulse();
    step_body_heat = body_->TakeStepHeat();
    const Point impulse = step_body_impulse / density;
    body_impulse_ += impulse;
    if (dt > 0) force = impulse / dt;
    if (t > 0) mean_force = body_impulse_ / t;
  }
  return {t, dt, t >= t_end, inner_particles, force, mean_force};
}

void CellCluster::InvalidateBodyCheck() {
  body_check_cells_.clear();
  body_check_bounds_.reset();
  body_check_cell_count_ = 0;
}

std::size_t CellCluster::CountBodyParticles() {
  if (!body_ || body_->PolygonCount() == 0) {
    InvalidateBodyCheck();
    return 0;
  }
  const auto bounds = body_->Bounds();
  const auto same = [](Point a, Point b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  };
  if (!body_check_bounds_ || body_check_cell_count_ != cells.size() ||
      !same(body_check_bounds_->first, bounds.first) ||
      !same(body_check_bounds_->second, bounds.second)) {
    // Publish a complete cache only after a successful rebuild.
    InvalidateBodyCheck();
    const auto [lower, upper] = bounds;
    for (auto& cell : cells) {
      const auto origin = cell.GetApex();
      const auto extent = cell.GetSize();
      if (origin.x <= upper.x && origin.x + extent.x >= lower.x &&
          origin.y <= upper.y && origin.y + extent.y >= lower.y &&
          origin.z <= upper.z && origin.z + extent.z >= lower.z)
        body_check_cells_.push_back(&cell);
    }
    body_check_cell_count_ = cells.size();
    body_check_bounds_ = bounds;
  }
  constexpr std::size_t chunk_size = 32;
  const auto chunks = body_check_cells_.size() / chunk_size +
                      (body_check_cells_.size() % chunk_size != 0);
  const auto workers = std::min(thread_pool_.GetThreadCount(), chunks);
  if (workers <= 1) {
    std::size_t count = 0;
    for (const auto* cell : body_check_cells_) count += cell->CountInnerParticles(*body_);
    return count;
  }
  std::atomic<std::size_t> next{0};
  std::vector<std::size_t> counts(workers, 0);
  std::vector<std::future<void>> futures;
  futures.reserve(workers);
  std::exception_ptr failure;
  try {
    for (std::size_t worker = 0; worker < workers; ++worker) {
      futures.emplace_back(thread_pool_.Execute([&, worker] {
        std::size_t count = 0;
        for (;;) {
          const auto begin = next.fetch_add(chunk_size, std::memory_order_relaxed);
          if (begin >= body_check_cells_.size()) break;
          const auto end = std::min(begin + chunk_size, body_check_cells_.size());
          for (auto i = begin; i < end; ++i)
            count += body_check_cells_[i]->CountInnerParticles(*body_);
        }
        counts[worker] = count;
      }));
    }
  } catch (...) { failure = std::current_exception(); }
  // Tasks refer to local counters and cluster data: drain before unwinding.
  for (auto& future : futures) {
    try { future.get(); }
    catch (...) { if (!failure) failure = std::current_exception(); }
  }
  if (failure) std::rethrow_exception(failure);
  std::size_t count = 0;
  for (auto local : counts) count += local;
  return count;
}

CellCluster::OutgoingSources CellCluster::GatherParticles() {
  auto describe = [&](bool include_cells) {
    OutgoingSources sources;
    std::size_t total = 0;
    auto append = [&](std::vector<Particle>& buffer) {
      if (buffer.empty()) return;
      if (buffer.size() > std::numeric_limits<std::size_t>::max() - total)
        throw std::length_error("Outgoing particle count overflow");
      total += buffer.size();
      sources.push_back({&buffer, total});
    };
    append(partile_buffer);
    if (include_cells) for (auto& cell : cells) append(cell.GetBuffer());
    return sources;
  };
  constexpr size_t chunk_size = 32;
  const auto chunks = cells.size() / chunk_size + (cells.size() % chunk_size != 0);
  const auto workers = std::min(thread_pool_.GetThreadCount(), chunks);
  const bool has_neighbors = std::any_of(cells.begin(), cells.end(),
      [](const Cell& cell) { return !cell.neighbors.empty(); });
  if (has_neighbors) {
    // Neighbor arrivals can affect cells later in this same traversal.
    for (auto& cell : cells) {
      cell.Sort();
      cell.SortNeighbors();
      auto& cell_buffer = cell.GetBuffer();
      partile_buffer.insert(partile_buffer.end(),
                            std::make_move_iterator(cell_buffer.begin()),
                            std::make_move_iterator(cell_buffer.end()));
      cell_buffer.clear();
    }
    return describe(false);
  }
  if (workers <= 1) {
    for (auto& cell : cells) cell.Sort();
    return describe(true);
  }

  std::atomic<size_t> next{0};
  auto sort = [&] {
    std::exception_ptr failure;
    for (;;) {
      const auto begin = next.fetch_add(chunk_size, std::memory_order_relaxed);
      if (begin >= cells.size()) break;
      const auto end = std::min(begin + chunk_size, cells.size());
      for (auto i = begin; i < end; ++i) {
        try { cells[i].Sort(); }
        catch (...) { if (!failure) failure = std::current_exception(); }
      }
    }
    if (failure) std::rethrow_exception(failure);
  };
  std::vector<std::future<void>> futures;
  futures.reserve(workers);
  std::exception_ptr failure;
  try {
    for (size_t worker = 0; worker < workers; ++worker)
      futures.emplace_back(thread_pool_.Execute(sort));
  } catch (...) { failure = std::current_exception(); }
  // Tasks access cell-owned buffers: drain all submissions before unwinding.
  for (auto& future : futures) {
    try { future.get(); }
    catch (...) { if (!failure) failure = std::current_exception(); }
  }
  if (failure) std::rethrow_exception(failure);

  return describe(true);
}

void CellCluster::BoundaryCondition() {
  for (auto& boundary : boundary_cond_outer) {
    // Keep face order: later faces also process arrivals emitted by earlier ones.
    boundary->ProcessOutgoing(partile_buffer, dt);
    auto inflow = boundary->GenerateInflow(dt);
    for (auto& batch : inflow.local)
      batch.cell.get().AppendGeneratedParticles(std::move(batch.particles));
    partile_buffer.insert(partile_buffer.end(),
                          std::make_move_iterator(inflow.transit.begin()),
                          std::make_move_iterator(inflow.transit.end()));
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
  if (buffer.empty()) return;
  DistributeSources({{&buffer, buffer.size()}}, buffer);
}

CellCluster::RoutingStorage CellCluster::DistributeSources(const OutgoingSources& sources,
                                    std::vector<Particle>& output) {
  if (sources.empty()) return {};
  const auto total = sources.back().end;
  const auto thread_count = std::max<std::size_t>(1, thread_pool_.GetThreadCount());
  const auto chunk_size = std::max<std::size_t>(1024,
      total / thread_count + (total % thread_count != 0));
  const auto chunk_count = total / chunk_size + (total % chunk_size != 0);
  std::vector<std::vector<Particle>> leftovers(chunk_count);
  std::vector<std::future<void>> futures;
  futures.reserve(chunk_count);
  auto distribute = [&](std::size_t begin, std::size_t end, std::vector<Particle>& remaining) {
    std::vector<std::size_t> accepted;
    accepted.reserve(std::min<std::size_t>(1024, end - begin));
    auto source = std::upper_bound(sources.begin(), sources.end(), begin,
        [](std::size_t index, const OutgoingSource& range) { return index < range.end; });
    while (begin < end) {
      const auto source_begin = source->end - source->particles->size();
      const auto stop = std::min(end, source->end);
      std::span<Particle> particles(*source->particles);
      std::size_t destination = 0;
      auto flush = [&] {
        if (!accepted.empty()) {
          cell_lookup_[destination]->AppendRoutedParticles(particles, accepted);
          accepted.clear();
        }
      };
      for (auto index = begin - source_begin; index < stop - source_begin; ++index) {
        std::size_t next_destination = 0;
        if (FindCellIndex(particles[index].position, next_destination) &&
            cell_lookup_[next_destination] &&
            cell_lookup_[next_destination]->ContainsPosition(particles[index].position)) {
          if (next_destination != destination) flush();
          destination = next_destination;
          accepted.push_back(index);
          if (accepted.size() == 1024) flush();
        } else {
          // Allocate only for actual boundary leftovers, never for the whole input.
          remaining.push_back(std::move(particles[index]));
        }
      }
      flush();
      begin = stop;
      ++source;
    }
  };
  if (chunk_count == 1) {
    distribute(0, total, leftovers.front());
  } else {
    std::exception_ptr failure;
    try {
      for (std::size_t chunk = 0; chunk < chunk_count; ++chunk) {
        const auto begin = chunk * chunk_size;
        const auto end = begin + std::min(chunk_size, total - begin);
        futures.emplace_back(thread_pool_.Execute([&, begin, end, chunk] {
          distribute(begin, end, leftovers[chunk]);
        }));
      }
    } catch (...) { failure = std::current_exception(); }
    for (auto& future : futures) {
      try { future.get(); }
      catch (...) { if (!failure) failure = std::current_exception(); }
    }
    if (failure) std::rethrow_exception(failure);
  }
  RoutingStorage storage;
  for (const auto& remaining : leftovers) {
    storage.size += remaining.size();
    storage.capacity += remaining.capacity();
  }
  // Reserve before clearing input. All tasks have finished; output may alias a source.
  output.reserve(storage.size);
  for (const auto& source : sources) source.particles->clear();
  output.clear();
  for (auto& remaining : leftovers)
    output.insert(output.end(), std::make_move_iterator(remaining.begin()),
                  std::make_move_iterator(remaining.end()));
  return storage;
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

  boundary->BindCells(cells);

  boundary_cond_outer.push_back(std::move(boundary));
}

void CellCluster::SetEndTime(double t_end) {
  this->data_t = 0;
  this->data_dt = t_end / 100 + 0.000001;
  this->t_end = t_end;
}
void CellCluster::CalculateCellTimeSteps() {
  constexpr size_t chunk_size = 32;
  const auto chunks = cells.size() / chunk_size + (cells.size() % chunk_size != 0);
  const auto workers = std::min(thread_pool_.GetThreadCount(), chunks);
  if (workers <= 1) {
    for (auto& cell : cells) cell.CalculateDt();
    return;
  }

  std::atomic<size_t> next{0};
  auto calculate = [&] {
    std::exception_ptr failure;
    for (;;) {
      const auto begin = next.fetch_add(chunk_size, std::memory_order_relaxed);
      if (begin >= cells.size()) break;
      const auto end = std::min(begin + chunk_size, cells.size());
      for (auto i = begin; i < end; ++i) {
        try { cells[i].CalculateDt(); }
        catch (...) { if (!failure) failure = std::current_exception(); }
      }
    }
    if (failure) std::rethrow_exception(failure);
  };
  std::vector<std::future<void>> futures;
  futures.reserve(workers);
  std::exception_ptr failure;
  try {
    for (size_t worker = 0; worker < workers; ++worker)
      futures.emplace_back(thread_pool_.Execute(calculate));
  } catch (...) { failure = std::current_exception(); }
  // Every submitted task must finish before publishing dt or unwinding.
  for (auto& future : futures) {
    try { future.get(); }
    catch (...) { if (!failure) failure = std::current_exception(); }
  }
  if (failure) std::rethrow_exception(failure);
}

void CellCluster::CalculateDt() {
  CalculateCellTimeSteps();
  // Read stored local results in the original order; no extra per-cell buffer.
  dt = std::numeric_limits<double>::infinity();

  cell_iter = cells.begin();
  for (auto& cell : cells) {
    double dtt = cell.PeekDt();
    if (dtt < dt) dt = dtt;
  }
  for (const auto& boundary : boundary_cond_outer)
    dt = std::min(dt, boundary->TransportTimeLimit());
  transport_dt = dt * Cu;
  collision_dt = std::numeric_limits<double>::infinity();
  if (collision_courant > 0) {
    for (const auto& cell : cells) {
      const double rate = cell.resolution.frequency;
      if (rate > 0) collision_dt = std::min(collision_dt, collision_courant / rate);
    }
  }
  dt = std::min({transport_dt, collision_dt, std::max(0., t_end - t)});

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

namespace mc3d {
void CellCluster::ApplyInitialStep(const SimulationConfig& config) {
  if(!config.initial_step) return;
  for(auto& cell:cells) {
    const auto& s=cell.GetCenter().x<config.initial_step->x?config.initial_step->left:config.initial_step->right;
    const double expected=density*cell.GetVolume()*s.density_ratio;
    if(!std::isfinite(expected) || expected>std::numeric_limits<unsigned int>::max()) throw std::invalid_argument("Initial population overflow");
    auto& random=cell.Random();
    const size_t count=static_cast<size_t>(expected)+(random.Uniform()<expected-std::floor(expected));
    cell.particles.clear();
    std::normal_distribution<double> normal(0,std::sqrt(s.temperature));
    auto origin=cell.GetCenter()-cell.GetSize()*.5; auto size=cell.GetSize();
    for(size_t i=0;i<count;++i) {
      Particle p;
      p.position=origin+Point(random.Uniform()*size.x,random.Uniform()*size.y,random.Uniform()*size.z);
      p.velocity={s.velocity[0]+normal(random.Engine()),s.velocity[1]+normal(random.Engine()),s.velocity[2]+normal(random.Engine())};
      cell.particles.push_back(p);
    }
    cell.SetTemperature(s.temperature);cell.SetVelocity({s.velocity[0],s.velocity[1],s.velocity[2]});
    cell.RefreshMoments();cell.GetDt();
  }
}
}
