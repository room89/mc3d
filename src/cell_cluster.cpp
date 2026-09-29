#include "cell_cluster.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <future>
#include <iterator>
#include <sstream>
#include <string>
#include <utils/logger.hpp>
#include <utils/stopwatch.hpp>

using namespace std;
using namespace mc3d;

namespace mc3d {
namespace {

constexpr std::size_t kEstimatedRowSize = 90;
constexpr int kOutputPrecision = 8;
constexpr double kSnapshotTolerance = 1e-9;

template <typename Int>
void AppendIntegral(std::string& out, Int value) {
  char buffer[32];
  auto result = std::to_chars(std::begin(buffer), std::end(buffer), value);
  if (result.ec == std::errc()) {
    out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
  } else {
    out.append(std::to_string(value));
  }
}

void AppendDouble(std::string& out, double value) {
  char buffer[64];
  auto result = std::to_chars(std::begin(buffer), std::end(buffer), value,
                              std::chars_format::general, kOutputPrecision);
  if (result.ec == std::errc()) {
    out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
  } else {
    std::ostringstream fallback;
    fallback.setf(std::ios::scientific);
    fallback.precision(kOutputPrecision > 0 ? kOutputPrecision - 1 : 0);
    fallback << value;
    out.append(fallback.str());
  }
}

}  // namespace

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

CellCluster::~CellCluster() { cells.clear(); }

bool CellCluster::Initialize(unsigned int ncx, unsigned int ncy,
                             unsigned int ncz, double density, double Kn,
                             double Cu, std::unique_ptr<Geometry>&& body,
                             double S, double alpha, double T) {
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

        temp_cell.SetApex(a);
        temp_cell.SetSize(dx, dy, dz);

        temp_cell.Initialize(N_, body_);

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

  this->density = N;

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

void CellCluster::SetBinaryOutput(bool enabled) { binary_output_ = enabled; }

void CellCluster::SetSnapshotInterval(double interval) {
  if (interval > 0.0) {
    snapshot_interval_ = interval;
    next_snapshot_time_ = interval;
  } else {
    snapshot_interval_ = 0.0;
    next_snapshot_time_ = 0.0;
  }
}

bool CellCluster::UpdateSnapshotSchedule() {
  if (snapshot_interval_ <= 0.0) {
    return true;
  }

  if (next_snapshot_time_ <= 0.0) {
    next_snapshot_time_ = snapshot_interval_;
  }

  if (t + kSnapshotTolerance >= next_snapshot_time_) {
    while (t + kSnapshotTolerance >= next_snapshot_time_) {
      next_snapshot_time_ += snapshot_interval_;
    }
    return true;
  }

  return false;
}

bool CellCluster::TimeStep() {
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
    for (auto& fut : futures) {
      fut.get();
    }
  }

  {
    utils::Stopwatch sw("TimeStep::GatherParticles");
    for (auto& cell : cells) {
      cell.CalculateVelocity();
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

  const bool should_write_snapshot = UpdateSnapshotSchedule();

  if (should_write_snapshot) {
    std::string file_name = "data";
    std::ostringstream ost;
    ost << t << ".dat";
    file_name += ost.str();
    {
      utils::Stopwatch sw("TimeStep::WriteFile");
      WriteFile(file_name);
    }
  }

  if (t >= t_end)
    return true;
  else
    return false;
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

bool CellCluster::WriteFile(const std::string& file_name) {
  std::ofstream file(file_name, std::ios::out | std::ios::binary);
  if (!file) {
    LOG_ERROR() << "Failed to open output file: " << file_name;
    return false;
  }

  if (binary_output_) {
    return WriteBinarySnapshot(file);
  }

  return WriteTextSnapshot(file);
}

bool CellCluster::WriteTextSnapshot(std::ofstream& file) {
  const std::size_t cell_count = cells.size();
  std::string buffer;
  buffer.reserve(1 + kEstimatedRowSize * cell_count);
  buffer.append("x;y;z;N;ro;T;vx;vy;vz;E\n");

  for (auto& cell : cells) {
    const auto center = cell.GetCenter();
    const auto velocity = cell.GetVelocity();
    const unsigned int particle_count = cell.GetParticleCount();
    const double density_value =
        cell.GetVolume() > 0 && density != 0.0
            ? static_cast<double>(particle_count) / (cell.GetVolume() * density)
            : 0.0;

    AppendDouble(buffer, center.x);
    buffer.push_back(';');
    AppendDouble(buffer, center.y);
    buffer.push_back(';');
    AppendDouble(buffer, center.z);
    buffer.push_back(';');
    AppendIntegral(buffer, particle_count);
    buffer.push_back(';');
    AppendDouble(buffer, density_value);
    buffer.push_back(';');
    AppendDouble(buffer, cell.GetTemperature());
    buffer.push_back(';');
    AppendDouble(buffer, velocity.x);
    buffer.push_back(';');
    AppendDouble(buffer, velocity.y);
    buffer.push_back(';');
    AppendDouble(buffer, velocity.z);
    buffer.push_back(';');
    AppendDouble(buffer, cell.GetEnergy());
    buffer.push_back('\n');
  }

  file.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  file.flush();
  return file.good();
}

bool CellCluster::WriteBinarySnapshot(std::ofstream& file) {
  static constexpr char kMagic[4] = {'M', 'C', '3', 'D'};
  file.write(kMagic, sizeof(kMagic));

  const std::uint32_t version = 1;
  file.write(reinterpret_cast<const char*>(&version), sizeof(version));

  const std::uint64_t record_count = static_cast<std::uint64_t>(cells.size());
  file.write(reinterpret_cast<const char*>(&record_count),
             sizeof(record_count));

  const double density_reference = density;
  file.write(reinterpret_cast<const char*>(&density_reference),
             sizeof(density_reference));

  if (!file) {
    return false;
  }

  for (auto& cell : cells) {
    const auto center = cell.GetCenter();
    const auto velocity = cell.GetVelocity();
    const std::uint32_t particle_count =
        static_cast<std::uint32_t>(cell.GetParticleCount());
    const double density_value =
        cell.GetVolume() > 0 && density != 0.0
            ? static_cast<double>(particle_count) / (cell.GetVolume() * density)
            : 0.0;
    const double temperature = cell.GetTemperature();
    const double energy = cell.GetEnergy();

    file.write(reinterpret_cast<const char*>(&center.x), sizeof(center.x));
    file.write(reinterpret_cast<const char*>(&center.y), sizeof(center.y));
    file.write(reinterpret_cast<const char*>(&center.z), sizeof(center.z));
    file.write(reinterpret_cast<const char*>(&particle_count),
               sizeof(particle_count));
    file.write(reinterpret_cast<const char*>(&density_value),
               sizeof(density_value));
    file.write(reinterpret_cast<const char*>(&temperature),
               sizeof(temperature));
    file.write(reinterpret_cast<const char*>(&velocity.x), sizeof(velocity.x));
    file.write(reinterpret_cast<const char*>(&velocity.y), sizeof(velocity.y));
    file.write(reinterpret_cast<const char*>(&velocity.z), sizeof(velocity.z));
    file.write(reinterpret_cast<const char*>(&energy), sizeof(energy));

    if (!file) {
      return false;
    }
  }

  file.flush();
  return file.good();
}

bool CellCluster::WriteSpeedFile(const char* file_name) {
  std::ofstream file1(file_name);
  cell_iter = cells.begin();
  Point ap;
  unsigned int N;
  for (auto& cell : cells) {
    ap = cell.GetCenter();
    N = cell.GetParticleCount();
    if (double(N) / (cell.GetVolume() * density) < 0) {
      bool a = body_->IsInnerPoint(ap);
      a = body_->IsInnerPoint(cell.GetApex());
    }
    file1 << ap << cell.GetVelocity() << '\n';
    cell_iter++;
  }
  file1.close();
  return true;
}

bool CellCluster::WriteFile() {
  std::ofstream file1("data.dat");
  cell_iter = cells.begin();
  unsigned int N;
  double av_den = 0;
  Point cell_size;
  double volume = 1;
  for (auto& cell : cells) {
    auto ap = cell.GetCenter();
    av_den += double(N) / np;
    file1 << ap << "\t"
          << double(cell.GetParticleCount()) / (cell.GetVolume() * density)
          << "\t" << cell.GetTemperature() << '\n';
    cell_iter++;
  }
  file1.close();
  return true;
}

bool CellCluster::WriteSpeedFile() {
  std::ofstream file1("speed.dat");
  cell_iter = cells.begin();
  Point ap;
  unsigned int N;
  double av_den = 0;
  for (auto& cell : cells) {
    ap = cell.GetCenter();
    N = cell.GetParticleCount();
    av_den += double(N) / np;
    file1 << ap << "\t" << cell.GetVelocity() << '\n';
    cell_iter++;
  }

  file1.close();
  return true;
}

bool CellCluster::WriteTimes() {
  if (proc_id == 0) {
    std::ofstream times_file("time.dat");

    deque<double>::iterator time_iter = times.begin();

    int i = 0;

    while (time_iter != times.end()) {
      times_file << i << "\t" << times[i] << "\t" << sort_times[i] << "\t"
                 << send_times[i] << "\t" << bound_times[i] << "\t"
                 << calc_times[i] << endl;

      i++;
      time_iter++;
    }

    times_file.close();
  }
  return 0;
}

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

void CellCluster::Compute() {
  size_t step = 0;
  while (t < t_end) {
    LOG_INFO() << "Time step:" << ++step;
    utils::Stopwatch stopwatch("ClusterStep");
    this->TimeStep();
  }
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

bool CellCluster::WriteCellFile(const std::string& init_file) {
  std::ofstream file(init_file);
  if (!file.is_open()) {
    LOG_ERROR() << "Can't open file: " << init_file;
    return false;
  }

  file << "cluster" << std::endl;

  cell_iter = cells.begin();
  for (auto& cell : cells) {
    Point V = cell.GetVelocity();
    double S = std::sqrt((V.x * V.x + V.y * V.y) / 2 * cell.GetTemperature());
    double alpha = M_PI / 2 - std::atan(V.x / V.y);
    if (V.x > 100000 * V.y) {
      V.x > 0 ? alpha = 0 : alpha = 3.141592654;
    }

    file << "\tcell" << endl;

    file << "\t\tapex\t" << cell.GetApex() << endl;

    file << "\t\tsize\t" << cell.GetSize() << endl;

    file << "\t\tS\t" << S << endl;

    file << "\t\talpha\t" << alpha << endl;

    file << "\t\tparticles\t" << cell.GetParticleCount() << endl;

    file << "\t\ttemperature\t" << cell.GetTemperature() << endl;

    file << "\t\tL\t" << cell.GetCharacteristicLength() << endl;

    file << "\tendcell" << endl;

    cell_iter++;
  }

  return true;
}

void CellCluster::SetDataSaveDtime(double data_dt) { this->data_dt = data_dt; }

void CellCluster::CleanInnerParticles() {
  for (auto& cell : cells) {
    cell.CleanInnerParticles(*body_);
  }
}
}  // namespace mc3d
