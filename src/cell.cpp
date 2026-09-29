#include "cell.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iterator>
#include <mutex>
#include <utils/logger.hpp>
#include <utils/utils.hpp>

namespace mc3d {
namespace {
const double Pi = 3.14159265358979;
}

Cell::Cell() : particles_mutex_(std::make_shared<std::mutex>()) {
  this->calc_time = 0;
  np = 1;
  body_mark = false;
  volume_ = 0.0;
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
  this->GenerateRandom(N);

  dt = 100000;
  CalculateTemperature();
  double c = sqrt(2 * T);
  double dtt = lx / (fabs(velocity.x) + c) < ly / (fabs(velocity.y) + c)
                   ? lx / (fabs(velocity.x) + c)
                   : ly / (fabs(velocity.y) + c);
  dtt = dtt < lz / (fabs(velocity.z) + c) ? dtt : lz / (fabs(velocity.z) + c);
  this->dt = dtt;
  auto& rng = utils::RandomEngine();
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

double Cell::GenerateRandom(size_t N) {
  if (N == 0) {
    return 0;
  }

  Point d(lx, ly, lz);

  std::vector<Particle> new_particles;
  new_particles.reserve(N);

  double ti = 0;
  const size_t nn2 = N / 2;

  for (size_t i = 0; i < nn2; ++i) {
    double rn1 = utils::Random01();
    double rn2 = utils::Random01();

    if (rn1 <= 0) rn1 = 1e-8;

    double slg = std::sqrt(2 * std::fabs(std::log(rn1)));

    Particle p1;
    Particle p2;

    const double cos_phi = std::cos(2 * Pi * rn2);
    const double sin_phi = std::sin(2 * Pi * rn2);

    p1.SetU(slg * cos_phi);
    p1.SetV(slg * sin_phi);

    p2.SetU(-slg * cos_phi);
    p2.SetV(-slg * sin_phi);

    double rn3 = utils::Random01();
    double rn4 = utils::Random01();

    if (rn3 <= 0) rn3 = 1e-8;

    slg = std::sqrt(2 * std::fabs(std::log(rn3)));

    const double cos_w = std::cos(2 * Pi * rn4);

    p1.SetW(slg * cos_w);
    p2.SetW(-slg * cos_w);

    ti += 2 * (p1.GetVelocity() * p1.GetVelocity());

    new_particles.push_back(p1);
    new_particles.push_back(p2);
  }

  if (N % 2 == 1) {
    Particle p;
    p.SetVelocity(Point(0, 0, 0));
    new_particles.push_back(p);
  }

  const double nt = 1.0 / static_cast<double>(N);
  ti = ti * nt / 3.0;
  if (ti <= 0) {
    ti = 1.0;
  }

  const double sf = std::sqrt(1.0 / ti);

  for (auto& particle : new_particles) {
    particle.velocity *= sf;
    particle.velocity = particle.velocity * std::sqrt(T) + velocity;

    const double rnx = utils::Random01();
    const double rny = utils::Random01();
    const double rnz = utils::Random01();

    const Point pos(d.x * rnx, d.y * rny, d.z * rnz);
    particle.position = apex + pos;
  }

  particles.reserve(particles.size() + new_particles.size());
  particles.insert(particles.end(), new_particles.begin(), new_particles.end());

  CalculateTemperature();

  return 0;
}

double Cell::GenerateRandom(unsigned int N, double T, Point V) {
  if (N == 0) {
    return 0;
  }

  Point d(lx, ly, lz);

  double ti = 0;

  std::vector<Particle> new_particles;
  new_particles.reserve(N);

  const unsigned int nn2 = N / 2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = utils::Random01();
    double rn2 = utils::Random01();

    if (rn1 <= 0) rn1 = 1e-5;

    double slg = std::sqrt(2 * std::fabs(std::log(rn1)));

    Particle p1;
    Particle p2;

    const double cos_phi = std::cos(2 * Pi * rn2);
    const double sin_phi = std::sin(2 * Pi * rn2);

    p1.SetU(slg * cos_phi);
    p1.SetV(slg * sin_phi);

    p2.SetU(-slg * cos_phi);
    p2.SetV(-slg * sin_phi);

    double rn3 = utils::Random01();
    double rn4 = utils::Random01();

    if (rn3 <= 0) rn3 = 1e-5;

    slg = std::sqrt(2 * std::fabs(std::log(rn3)));

    const double cos_w = std::cos(2 * Pi * rn4);

    p1.SetW(slg * cos_w);
    p2.SetW(-slg * cos_w);

    ti += 2 * (p1.GetVelocity() * p1.GetVelocity());

    new_particles.push_back(p1);
    new_particles.push_back(p2);
  }

  if (N % 2 == 1) {
    Particle p;
    p.SetVelocity(Point(0, 0, 0));
    new_particles.push_back(p);
  }

  const double nt = 1.0 / static_cast<double>(N);
  ti = ti * nt / 3.0;
  if (ti <= 0) ti = 1.0;

  const double sf = std::sqrt(1.0 / ti);
  const double sqrt_T = std::sqrt(T);

  for (auto& particle : new_particles) {
    particle.velocity *= sf;
    particle.velocity *= sqrt_T;
    particle.velocity += V;

    const double rnx = utils::Random01();
    const double rny = utils::Random01();
    const double rnz = utils::Random01();

    particle.position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);
  }

  particles.reserve(particles.size() + new_particles.size());
  particles.insert(particles.end(), new_particles.begin(), new_particles.end());

  auto& rng = utils::RandomEngine();
  std::shuffle(particles.begin(), particles.end(), rng);

  return 0;
}

double Cell::GenerateFreeRandom(unsigned int N, double T, Point V, Point nrml) {
  if (N == 0) {
    return 0;
  }

  Point d(lx, ly, lz);

  double ti = 0;
  std::vector<Particle> added_particles;
  added_particles.reserve(N);
  Point dvel(0, 0, 0);

  unsigned int i = 0;

  while (i < N) {
    double rn1 = utils::Random01();
    double rn2 = utils::Random01();

    if (rn1 <= 0) rn1 = 1e-5;

    double slg = std::sqrt(2 * T * std::fabs(std::log(rn1)));

    Particle p1;
    const double cos_phi = std::cos(2 * Pi * rn2);
    const double sin_phi = std::sin(2 * Pi * rn2);

    p1.SetU(slg * cos_phi);
    p1.SetV(slg * sin_phi);

    double rn3 = utils::Random01();
    double rn4 = utils::Random01();

    if (rn3 <= 0) rn3 = 1e-5;

    slg = std::sqrt(2 * T * std::fabs(std::log(rn3)));
    const double cos_w = std::cos(2 * Pi * rn4);

    p1.SetW(slg * cos_w);

    if ((p1.GetVelocity() + V) * nrml < 0) {
      added_particles.push_back(p1);
      ti += (p1.velocity * p1.velocity);
      dvel += p1.velocity;
      ++i;
    }
  }

  const double nt = 1.0 / static_cast<double>(N);
  dvel *= nt;
  ti = (ti * nt - dvel * dvel) / 3;
  if (ti <= 0) ti = T;

  const double sf = std::sqrt(T / ti);

  for (auto& particle : added_particles) {
    particle.velocity = (particle.velocity - dvel) * sf + V;
    if (particle.velocity * nrml < 0) {
      particle.velocity = 2 * V - particle.velocity;
    }
  }

  for (auto& particle : added_particles) {
    const double rnx = utils::Random01();
    const double rny = utils::Random01();
    const double rnz = utils::Random01();

    particle.position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);
  }

  particles.reserve(particles.size() + added_particles.size());
  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

double Cell::GenerateHyperFreeRandom(unsigned int N, Point V, double T) {
  std::vector<Particle> added_particles;
  Point position;
  Point d(lx, ly, lz);
  double A = sqrt(3 * T);

  added_particles.reserve(N);

  for (unsigned int i = 0; i < N; i++) {
    double rnx = utils::RandomDouble(0, 1);
    double rny = utils::RandomDouble(0, 1);
    double rnz = utils::RandomDouble(0, 1);

    position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);

    Particle new_particle;

    new_particle.position = position;
    Point noise(A * utils::RandomDouble(-1, 1), A * utils::RandomDouble(-1, 1),
                A * utils::RandomDouble(-1, 1));
    new_particle.velocity = V + noise;

    added_particles.push_back(new_particle);
  }

  Point av_velocity;
  for (const auto& particle : added_particles) {
    av_velocity += particle.velocity;
  }
  av_velocity /= particles.size();

  particles.reserve(particles.size() + added_particles.size());
  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

void Cell::Collisions() {
  if (particles.size() <= 5) return;

  CalculateKn();

  double g_max = 2 * sqrt(CalculateTemperature());
  double factor = 2 * sqrt(2.) * L * Kn_l / particles.size();
  double frequency_t = factor / g_max;

  double t = 0;

  auto& rng = utils::RandomEngine();

  auto particle_1 = particles.begin();
  auto particle_2 = particles.begin();

  particle_2++;

  while (t <= dt) {
    double r = utils::Random01();

    if (r <= 0) r = 3.e-5;

    double tau = -frequency_t * log(r);
    t += tau;
    if (t > dt) break;

    Point v1 = particle_1->GetVelocity();
    Point v2 = particle_2->GetVelocity();

    Point velocity_sum = v1 + v2;

    Point g = v2 - v1;
    double gmod = g.Mod();

    if (g_max < gmod) {
      g_max = gmod;
      t -= tau;
      frequency_t = factor / g_max;
      continue;
    }

    double rr = utils::Random01();

    if (gmod / g_max > rr)  // g == gmod?
    {
      double r1 = utils::Random01();
      double r2 = utils::Random01();

      double Fi = 2. * Pi * r1;
      double Eta = 2. * acos(r2);

      double gx = g.x;
      double gy = g.y;
      double gz = g.z;
      double gxz = sqrt(gx * gx + gz * gz);

      Point g1;

      //                              Point g1(gmod * sin(Pi * r1) * cos(2. * Pi
      //                              * r2),
      //                                       gmod * sin(Pi * r1) * sin(2. * Pi
      //                                       * r2), gmod * cos(Pi * r1));

      if (gxz > 1.e-6) {
        g1.x = gx * cos(Eta) -
               sin(Eta) * (gmod * gz * cos(Fi) - gx * gy * sin(Fi)) / gxz;
        g1.y = gy * cos(Eta) - gxz * sin(Fi) * sin(Eta);
        g1.z = gz * cos(Eta) +
               sin(Eta) * (gmod * gx * cos(Fi) + gy * gz * sin(Fi)) / gxz;
      } else {
        g1.x = -sin(Eta) * gmod * (cos(Fi) - sin(Fi)) / std::sqrt(2.);
        g1.y = gmod * cos(Eta);
        g1.z = sin(Eta) * gmod * (cos(Fi) + sin(Fi)) / std::sqrt(2.);
      }

      v1 = (velocity_sum - g1) * 0.5;
      v2 = (velocity_sum + g1) * 0.5;

      particle_1->velocity = v1;
      particle_2->velocity = v2;
    }

    particle_2++;

    if (particle_2 == particles.end()) {
      std::shuffle(particles.begin(), particles.end(), rng);

      particle_1 = particles.begin();
      particle_2 = particles.begin();
      particle_2++;
    } else {
      particle_1 = particle_2;
      particle_2++;
      if (particle_2 == particles.end()) {
        std::shuffle(particles.begin(), particles.end(), rng);

        particle_1 = particles.begin();
        particle_2 = particles.begin();
        particle_2++;
      }
    }
  }
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

bool Cell::TryAcceptParticle(Particle& particle) {
  const bool inside_x =
      (particle.position.x >= apex.x) && (particle.position.x < apex.x + lx);
  const bool inside_y =
      (particle.position.y >= apex.y) && (particle.position.y < apex.y + ly);
  const bool inside_z =
      (particle.position.z >= apex.z) && (particle.position.z < apex.z + lz);

  if (!(inside_x && inside_y && inside_z)) {
    return false;
  }

  auto mutex_ptr = particles_mutex_;
  if (!mutex_ptr) {
    mutex_ptr = std::make_shared<std::mutex>();
    particles_mutex_ = mutex_ptr;
  }

  std::lock_guard<std::mutex> lock(*mutex_ptr);
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

double Cell::GetTemperatureRaw() { return T; }

void Cell::SetParameters(double S, double alpha, double T) {
  this->T = T;
  velocity =
      Point(S * sqrt(2 * T) * cos(alpha), S * sqrt(2 * T) * sin(alpha), 0);
}

void Cell::SetVelocity(Point velocity) { this->velocity = velocity; }

void Cell::SetTemperature(double t) { this->T = t; }

double Cell::GetKn() { return Kn; }

double Cell::CalculateKn() {
  Kn_l = Kn * np / particles.size();
  return Kn;
}

double Cell::GetTemperature() {
  CalculateTemperature();
  return T;
}

double Cell::GetU() {
  velocity = Point(0, 0, 0);
  for (const auto& particle : particles) {
    velocity += particle.velocity;
  }
  if (!particles.empty()) {
    velocity /= static_cast<double>(particles.size());
  }
  return velocity.x;
}

double Cell::GetV() {
  velocity = Point(0, 0, 0);
  for (const auto& particle : particles) {
    velocity += particle.velocity;
  }
  if (!particles.empty()) {
    velocity /= static_cast<double>(particles.size());
  }
  return velocity.y;
}

double Cell::GetW() {
  velocity = Point(0, 0, 0);
  for (const auto& particle : particles) {
    velocity += particle.velocity;
  }
  if (!particles.empty()) {
    velocity /= static_cast<double>(particles.size());
  }
  return velocity.z;
}

double Cell::GetEnergy() {
  double E = 0;
  for (const auto& particle : particles) {
    E += particle.velocity * particle.velocity;
  }
  if (!particles.empty()) {
    E /= 2 * particles.size();
  }
  return E;
}

Point Cell::CalculateVelocity() const {
  auto velocity = Point(0, 0, 0);
  for (const auto& particle : particles) {
    velocity += particle.GetVelocity();
  }

  if (!particles.empty()) {
    velocity /= static_cast<double>(particles.size());
  }

  return velocity;
}

double Cell::CalculateTemperature() {
  if (particles.empty()) {
    velocity = Point(0, 0, 0);
    T = 0;
    return T;
  }

  Point av_vel(0, 0, 0);
  double E = 0;

  for (const auto& particle : particles) {
    av_vel += particle.velocity;
    E += particle.velocity * particle.velocity;
  }

  const double inv_count = 1.0 / static_cast<double>(particles.size());
  av_vel *= inv_count;
  E *= inv_count;
  velocity = av_vel;
  T = (E - velocity * velocity) / 3;
  return T;
}

double Cell::CalculateDt() {
  CalculateTemperature();
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
  if (body_mark)
    body_boundary.BoundaryCondition(particles, dt);
  else
    MoveParticles();
}

void Cell::AttachThreadMark(int& ptr) { thread_mark = ptr; }

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
  volume_ = body_boundary.CalcCellVolume(apex, GetSize(), &mass_center);
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

Point Cell::GetVelocity() {
  velocity = Point(0, 0, 0);
  for (const auto& particle : particles) {
    velocity += particle.velocity;
  }
  if (!particles.empty()) {
    velocity /= static_cast<double>(particles.size());
  }
  return velocity;
}

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

bool operator>(const Cell& a, const Cell& b) {
  return a.calc_time > b.calc_time;
}

bool operator<(const Cell& a, const Cell& b) {
  return b.calc_time > a.calc_time;
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