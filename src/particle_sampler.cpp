#include "particle_sampler.h"
#include <algorithm>
#include <cmath>
#include <limits>
namespace mc3d {
namespace { constexpr double Pi = 3.14159265358979; }
double ParticleSampler::Initial(std::vector<Particle>& particles, Point apex, Point size, size_t N, double T, Point velocity, utils::RandomStream& random) {
  if (N == 0) {
    return 0;
  }

  Point d = size;

  std::vector<Particle> new_particles;
  new_particles.reserve(N);

  double ti = 0;
  const size_t nn2 = N / 2;

  for (size_t i = 0; i < nn2; ++i) {
    double rn1 = random.Uniform();
    double rn2 = random.Uniform();

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

    double rn3 = random.Uniform();
    double rn4 = random.Uniform();

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

    const double rnx = random.Uniform();
    const double rny = random.Uniform();
    const double rnz = random.Uniform();

    const Point pos(d.x * rnx, d.y * rny, d.z * rnz);
    particle.position = apex + pos;
  }

  particles.reserve(particles.size() + new_particles.size());
  particles.insert(particles.end(), new_particles.begin(), new_particles.end());


  return 0;
}

double ParticleSampler::Reservoir(std::vector<Particle>& particles, Point apex, Point size, unsigned int N, double T, Point V, utils::RandomStream& random) {
  if (N == 0) {
    return 0;
  }

  Point d = size;

  double ti = 0;

  std::vector<Particle> new_particles;
  new_particles.reserve(N);

  const unsigned int nn2 = N / 2;

  for (unsigned int i = 0; i < nn2; i++) {
    double rn1 = random.Uniform();
    double rn2 = random.Uniform();

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

    double rn3 = random.Uniform();
    double rn4 = random.Uniform();

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

    const double rnx = random.Uniform();
    const double rny = random.Uniform();
    const double rnz = random.Uniform();

    particle.position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);
  }

  particles.reserve(particles.size() + new_particles.size());
  particles.insert(particles.end(), new_particles.begin(), new_particles.end());

  auto& rng = random.Engine();
  std::shuffle(particles.begin(), particles.end(), rng);

  return 0;
}

double ParticleSampler::Free(std::vector<Particle>& particles, Point apex, Point size, unsigned int N, double T, Point V, Point nrml, utils::RandomStream& random) {
  if (N == 0) {
    return 0;
  }

  Point d = size;

  double ti = 0;
  std::vector<Particle> added_particles;
  added_particles.reserve(N);
  Point dvel(0, 0, 0);

  unsigned int i = 0;

  while (i < N) {
    double rn1 = random.Uniform();
    double rn2 = random.Uniform();

    if (rn1 <= 0) rn1 = 1e-5;

    double slg = std::sqrt(2 * T * std::fabs(std::log(rn1)));

    Particle p1;
    const double cos_phi = std::cos(2 * Pi * rn2);
    const double sin_phi = std::sin(2 * Pi * rn2);

    p1.SetU(slg * cos_phi);
    p1.SetV(slg * sin_phi);

    double rn3 = random.Uniform();
    double rn4 = random.Uniform();

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
    const double rnx = random.Uniform();
    const double rny = random.Uniform();
    const double rnz = random.Uniform();

    particle.position = apex + Point(d.x * rnx, d.y * rny, d.z * rnz);
  }

  particles.reserve(particles.size() + added_particles.size());
  particles.insert(particles.end(), added_particles.begin(),
                   added_particles.end());

  return 0;
}

void ParticleSampler::HyperFree(Point apex, Point size, unsigned int N, Point V, double T, Point normal, double dt, InnerBoundary& body_boundary, std::vector<Particle>& arrivals, utils::RandomStream& random) {
  std::normal_distribution<double> thermal(0, std::sqrt(T));
  const double inward_drift = -(V * normal);
  const Point tangent_drift = V - normal * (V * normal);
  std::vector<Particle> particle(1);
  for (unsigned int i = 0; i < N; ++i) {
    Point noise(thermal(random.Engine()), thermal(random.Engine()),
                thermal(random.Engine()));
    const double w = utils::IncomingNormalSpeed(inward_drift, T, random);
    particle[0].velocity = tangent_drift + noise - normal * (noise * normal) - normal * w;
    Point position = apex + Point(size.x * random.Uniform(), size.y * random.Uniform(),
                                  size.z * random.Uniform());
    auto face = [](double a, double size, double n) {
      return std::nextafter(n < 0 ? a : a + size,
          n < 0 ? std::numeric_limits<double>::infinity()
                : -std::numeric_limits<double>::infinity());
    };
    if (normal.x != 0) position.x = face(apex.x, size.x, normal.x);
    if (normal.y != 0) position.y = face(apex.y, size.y, normal.y);
    if (normal.z != 0) position.z = face(apex.z, size.z, normal.z);
    particle[0].position = position;
    // Uniform entry time within this step; advance only the remaining time.
    const double remaining = dt * random.Uniform();
    if (body_boundary.GetGeometryPtr() && body_boundary.GetGeometryPtr()->PolygonCount() > 0)
      body_boundary.BoundaryCondition(particle, remaining, random);
    else
      particle[0].position += particle[0].velocity * remaining;
    // The arrival may cross into another cell (or leave at a corner).
    arrivals.push_back(particle[0]);
  }
}

}
