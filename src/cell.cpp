#include "cell.h"
#include "particle_sampler.h"
#include "collision_model.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <utils/logger.hpp>
#include <utils/utils.hpp>

namespace mc3d {

Cell::Cell(std::shared_ptr<utils::RandomStream> random) {
  if (random) random_ = std::move(random);
  body_mark = false;
  volume_ = 0.0;
}

Cell::Cell(uint32_t seed, uint64_t stream_id) : Cell() {
  random_.emplace<utils::RandomStream>(seed, stream_id);
}

Cell::~Cell() {
  particles.clear();
  particle_buffer.clear();
  neighbors.clear();
}

void Cell::SetSize(double lx, double ly, double lz) {
  if ((lx <= 0) || (ly <= 0) || (lz <= 0)) return;
  this->lx = lx;
  this->ly = ly;
  this->lz = lz;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
  volume_ = std::abs(lx * ly * lz);
}

void Cell::SetSize(Point dl) {
  if ((dl.x <= 0) || (dl.y <= 0) || (dl.z <= 0)) return;
  lx = dl.x;
  ly = dl.y;
  lz = dl.z;
  L = lx > ly ? lx : ly;
  L = L > lz ? L : lz;
  volume_ = std::abs(lx * ly * lz);
}

void Cell::SetApex(Point a) { apex = a; }

bool Cell::Initialize(size_t N, const std::unique_ptr<Geometry>& bbody) {
  // Keep the reservoir normalization, not the fluctuating current population.
  // Empty standalone cells retain the explicit/default normalization.
  if (N > 0) SetReferenceParticleCount(static_cast<double>(N));
  this->GenerateRandom(N);

  dt = 100000;
  RefreshMoments();
  double c = sqrt(2 * T);
  double dtt = lx / (fabs(velocity.x) + c) < ly / (fabs(velocity.y) + c)
                   ? lx / (fabs(velocity.x) + c)
                   : ly / (fabs(velocity.y) + c);
  dtt = dtt < lz / (fabs(velocity.z) + c) ? dtt : lz / (fabs(velocity.z) + c);
  this->dt = dtt;
  auto& rng = Random().Engine();
  std::shuffle(particles.begin(), particles.end(), rng);

  if (bbody) {
    body_boundary.SetGeometry(*bbody);
    body_mark = body_boundary.AddPolygon(*bbody, GetCenter(), L);
    if (body_mark) {
      CalculateVolume();
    }
  } else {
    body_mark = false;
  }

  return true;
}

unsigned int Cell::GetParticleCount() const { return particles.size(); }

std::size_t Cell::CountInnerParticles(const Geometry& body) const {
  return std::count_if(particles.begin(), particles.end(),
                       [&body](const Particle& particle) {
                         return body.IsInnerPoint(particle.position);
                       });
}

double Cell::GenerateRandom(size_t N) {
  const auto result = ParticleSampler::Initial(particles, apex, GetSize(), N, T, velocity, Random());
  if (N > 0) RefreshMoments();
  return result;
}
double Cell::GenerateRandom(unsigned int N, double temperature, Point drift) {
  return ParticleSampler::Reservoir(particles, apex, GetSize(), N, temperature, drift, Random());
}
std::vector<Particle> Cell::SampleFreeInflow(unsigned int N, double temperature, Point drift, Point normal) {
  std::vector<Particle> arrivals;
  ParticleSampler::Free(arrivals, apex, GetSize(), N, temperature, drift, normal, Random());
  return arrivals;
}
void Cell::AppendGeneratedParticles(std::vector<Particle>&& incoming) {
  particles.insert(particles.end(), std::make_move_iterator(incoming.begin()),
                   std::make_move_iterator(incoming.end()));
}
void Cell::GenerateHyperFreeRandom(unsigned int N, Point drift, double temperature,
                                  Point normal, double step, std::vector<Particle>& arrivals) {
  ParticleSampler::HyperFree(apex, GetSize(), N, drift, temperature, normal, step,
                            body_boundary, arrivals, Random());
}

void Cell::Collisions() {
  if (particles.size() <= 5) return;
  RefreshMoments();
  CollisionsWithPreparedMoments();
}

void Cell::CollisionsWithPreparedMoments() {
  if (particles.size() <= 5) return;
  CalculateKn();
  CollisionModel::Apply(particles, velocity, L, Kn_l, dt, Random());
}

void Cell::MoveParticles() {
  for (auto& particle : particles) {
    particle.position += particle.velocity * dt;
  }
}

void Cell::Sort() {
  auto inside = [this](const Particle& particle) {
    const bool inside_x =
        (particle.position.x >= apex.x) && (particle.position.x < apex.x + lx);
    const bool inside_y =
        (particle.position.y >= apex.y) && (particle.position.y < apex.y + ly);
    const bool inside_z =
        (particle.position.z >= apex.z) && (particle.position.z < apex.z + lz);
    return inside_x && inside_y && inside_z;
  };

  auto first_outside =
      std::partition(particles.begin(), particles.end(), inside);

  for (auto it = first_outside; it != particles.end(); ++it) {
    particle_buffer.push_back(std::move(*it));
  }

  particles.erase(first_outside, particles.end());
}

void Cell::SetCharacteristicLength(double L) { this->L = L; }

void Cell::AddParticle(std::vector<Particle>& incoming) {
  auto write_it = incoming.begin();

  for (auto it = incoming.begin(); it != incoming.end(); ++it) {
    const bool inside_x =
        (it->position.x >= apex.x) && (it->position.x < apex.x + lx);
    const bool inside_y =
        (it->position.y >= apex.y) && (it->position.y < apex.y + ly);
    const bool inside_z =
        (it->position.z >= apex.z) && (it->position.z < apex.z + lz);

    if (inside_x && inside_y && inside_z) {
      particles.push_back(std::move(*it));
    } else {
      if (write_it != it) {
        *write_it = std::move(*it);
      }
      ++write_it;
    }
  }

  incoming.erase(write_it, incoming.end());
}

void Cell::AppendRoutedParticles(std::span<Particle> source,
                                 std::span<const std::size_t> indices) {
  if (indices.empty()) return;
  std::lock_guard<std::mutex> lock(particles_mutex_);
  if (indices.size() > particles.max_size() - particles.size())
    throw std::length_error("Particle batch exceeds cell capacity");
  const auto required = particles.size() + indices.size();
  if (required > particles.capacity()) {
    // Keep amortized growth: exact reserves would reallocate on every small influx.
    const auto capacity = particles.capacity();
    const auto grown = capacity + std::min(capacity, particles.max_size() - capacity);
    particles.reserve(std::max(required, grown));
  }
  for (auto index : indices) particles.push_back(std::move(source[index]));
}

bool Cell::TryAcceptParticle(Particle& particle) {
  if (!ContainsPosition(particle.position)) return false;

  std::lock_guard<std::mutex> lock(particles_mutex_);
  particles.push_back(std::move(particle));
  return true;
}

std::vector<Particle>& Cell::GetBuffer() { return particle_buffer; }

/*var Cell::get_var()
{
        var a(velocity.x, velocity.y, velocity.z, T, n);
        return a;
}*/

Point Cell::GetApex() const { return apex; }

Point Cell::GetCenter() const { return apex + this->GetSize() / 2; }

Point Cell::GetMassCenter() const { return mass_center; }

double Cell::GetDt() {
  this->CalculateDt();
  return dt;
}

double Cell::PeekDt() const { return dt; }

void Cell::SetDt(double dt) { this->dt = dt; }

void Cell::SetParameters(double S, double alpha, double T) {
  this->T = T;
  velocity =
      Point(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
}

void Cell::SetVelocity(Point velocity) { this->velocity = velocity; }

void Cell::SetTemperature(double t) { this->T = t; }

double Cell::GetKn() { return Kn; }

void Cell::SetReferenceParticleCount(double count) {
  if (!(count > 0) || !std::isfinite(count))
    throw std::invalid_argument("Reference particle count must be finite and positive");
  reference_particle_count_ = count;
}

double Cell::CalculateKn() {
  if (particles.empty()) return Kn_l = std::numeric_limits<double>::infinity();
  // n/n0 = N / (N0 * Vgas/Vcell); N0 includes the model-particle weight.
  Kn_l = Kn * reference_particle_count_ * GetVolume() /
         (lx * ly * lz * particles.size());
  return Kn_l;
}

CellMoments Cell::ComputeMoments() const {
  CellMoments result;
  result.count = particles.size();
  double squared_speed = 0;
  for (const auto& particle : particles) {
    result.velocity += particle.velocity;
    squared_speed += particle.velocity * particle.velocity;
  }
  if (result.count != 0) {
    const double inv_count = 1.0 / static_cast<double>(result.count);
    result.velocity *= inv_count;
    squared_speed *= inv_count;
    result.energy = squared_speed / 2;
    result.temperature = (squared_speed - result.velocity * result.velocity) / 3;
  }
  return result;
}

void Cell::RefreshMoments() {
  const auto moments = ComputeMoments();
  velocity = moments.velocity;
  T = moments.temperature;
}

double Cell::CalculateDt() {
  RefreshMoments();
  dt = 1000000;
  double c = sqrt(2 * T);
  double dtt = min(lx / (fabs(velocity.x) + c), ly / (fabs(velocity.y) + c));
  dtt = min(dtt, lz / (fabs(velocity.z) + c));
  if (dtt < dt) dt = dtt;
  // if(dt < 0.00001) cin>>c;
  return dt;
}

void Cell::SetKn(double Kn) { this->Kn = Kn; }

void Cell::AddNeighbor(Cell& neighbor) { neighbors.emplace_back(neighbor); }
void Cell::SortNeighbors() {
  for (Cell& neighbor : neighbors) {
    neighbor.AddParticle(particle_buffer);
  }
}

void Cell::Calculate() {
  Collisions();
  AdvanceParticles();
}

void Cell::CalculateWithPreparedMoments() {
  CollisionsWithPreparedMoments();
  AdvanceParticles();
}

void Cell::AdvanceParticles() {
  if (body_boundary.GetGeometryPtr() &&
      body_boundary.GetGeometryPtr()->PolygonCount() > 0)
    body_boundary.BoundaryCondition(particles, dt, Random());
  else
    MoveParticles();
}


Point Cell::GetSize() const { return Point(abs(lx), abs(ly), abs(lz)); }

double Cell::GetCharacteristicLength() { return L; }

std::deque<std::unique_ptr<Cell>> Cell::Fragment(
    const std::unique_ptr<Geometry>& body) {
  std::deque<std::unique_ptr<Cell>> new_cells;
  for (size_t i = 0; i < 8; i++) {
    new_cells.emplace_back(std::make_unique<Cell>());
  }

  Point a(0, 0, 0);
  Point dl(lx / 2, ly / 2, lz / 2);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->SetSize(dl);
  }

  a = apex;
  new_cells[0]->SetApex(a);

  Point shift(lx / 2, 0, 0);

  a = apex + shift;
  new_cells[1]->SetApex(a);

  shift.Set(0, ly / 2, 0);
  a = apex + shift;
  new_cells[2]->SetApex(a);

  shift.Set(0, 0, lz / 2);
  a = apex + shift;
  new_cells[3]->SetApex(a);

  shift.Set(lx / 2, ly / 2, 0);
  a = apex + shift;
  new_cells[4]->SetApex(a);

  shift.Set(lx / 2, 0, lz / 2);
  a = apex + shift;
  new_cells[5]->SetApex(a);

  shift.Set(0, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[6]->SetApex(a);

  shift.Set(lx / 2, ly / 2, lz / 2);
  a = apex + shift;
  new_cells[7]->SetApex(a);

  for (size_t i = 0; i < 8; i++) {
    new_cells[i]->Initialize(particles.size() / 8, body);
  }

  return new_cells;
}

void Cell::SetInnerBoundary(InnerBoundary bound) { body_boundary = bound; }

double Cell::CalculateVolume() {
  volume_ = body_boundary.CalcCellVolume(apex, GetSize(), Random(), &mass_center);
  if (volume_ < 0) {
    return volume_ = GetSize().Volume();
  }
  return volume_;
}

double Cell::GetVolume() const {
  if (volume_ > 0.0) {
    return volume_;
  }
  return std::abs(lx * ly * lz);
}

bool Cell::GetBodyMark() { return body_mark; }

Point Cell::GetParticleMassCenter() {
  Point particle_mass_center(0, 0, 0);
  for (const auto& particle : particles) {
    particle_mass_center += particle.position;
  }
  if (!particles.empty()) {
    particle_mass_center /= static_cast<double>(particles.size());
  }
  return particle_mass_center;
}

void Cell::CleanInnerParticles(const Geometry& body) {
  if (particles.size() < 100) {
    Point ad = Point(0, 0, 0);
  }
  particle_buffer.insert(particle_buffer.end(),
                         std::make_move_iterator(particles.begin()),
                         std::make_move_iterator(particles.end()));
  const size_t buffered_count = particle_buffer.size();
  size_t removed_count = 0;
  std::array<Point, 3> sample_positions;
  std::array<Point, 3> sample_velocities;
  size_t sample_count = 0;

  particles.clear();

  for (auto& particle : particle_buffer) {
    if (!body.IsInnerPoint(particle.position)) {
      particles.push_back(std::move(particle));
    } else {
      if (sample_count < sample_positions.size()) {
        sample_positions[sample_count] = particle.position;
        sample_velocities[sample_count] = particle.velocity;
        ++sample_count;
      }
      ++removed_count;
    }
  }
  if (removed_count > 0) {
    LOG_WARNING() << "Removed " << removed_count
                  << " particles that entered the body geometry. Cell center: "
                  << GetCenter() << " size: " << GetSize()
                  << " buffer size: " << buffered_count;
    for (size_t idx = 0; idx < sample_count; ++idx) {
      LOG_DEBUG() << "Removed particle #" << idx
                  << " position=" << sample_positions[idx]
                  << " velocity=" << sample_velocities[idx];
    }
  }
  particle_buffer.clear();
}

void Cell::SetBodyMark(bool mark) { this->body_mark = mark; }

bool Cell::DebugTestParticle() {
  int i = 0;

  auto* geometry = this->body_boundary.GetGeometryPtr();
  if (!geometry) {
    return false;
  }

  for (const auto& particle : particles) {
    if (geometry->IsInnerPoint(particle.GetPosition())) {
      LOG_DEBUG() << particle.GetPosition() << "\t" << particle.GetVelocity();
      geometry->IsInnerPoint(particle.GetPosition());
      i++;
    }
  }
  return i > 0;
}
}  // namespace mc3d
