#include "collision_model.h"
#include <algorithm>
#include <cmath>
namespace mc3d {
namespace { constexpr double Pi = 3.14159265358979; }
void CollisionModel::Apply(std::vector<Particle>& particles, Point velocity, double L, double Kn_l, double dt, utils::RandomStream& random) {
  if (particles.size() <= 5) return;


  // Bound every relative speed by twice the largest peculiar speed. A thermal
  // mean is not an upper bound and causes population-dependent missed events.
  double max_speed_squared = 0;
  auto speed_squared = [velocity](const Point& v) {
    const Point c = v - velocity;
    return c * c;
  };
  for (const auto& p : particles)
    max_speed_squared = std::max(max_speed_squared, speed_squared(p.velocity));
  if (max_speed_squared == 0) return;
  double g_max = 2 * std::sqrt(max_speed_squared);
  // N*(N-1)/2 distinct unordered pairs, with particle weight set by N0.
  double factor = 2 * sqrt(2.) * L * Kn_l / (particles.size() - 1);
  double frequency_t = factor / g_max;

  double t = 0;

  auto& rng = random.Engine();

  std::uniform_int_distribution<std::size_t> first(0, particles.size()-1);
  std::uniform_int_distribution<std::size_t> second(0, particles.size()-2);

  while (t <= dt) {
    double r = random.Uniform();

    if (r <= 0) r = 3.e-5;

    double tau = -frequency_t * log(r);
    t += tau;
    if (t > dt) break;

    // Every attempt must sample all particles, even when a time step contains
    // fewer attempts than N/2. Restarting at the vector prefix biases collisions.
    const auto i = first(rng);
    auto j = second(rng);
    if (j >= i) ++j;
    auto particle_1 = particles.begin() + i;
    auto particle_2 = particles.begin() + j;

    Point v1 = particle_1->GetVelocity();
    Point v2 = particle_2->GetVelocity();

    Point velocity_sum = v1 + v2;

    Point g = v2 - v1;
    double gmod = g.Mod();

    double rr = random.Uniform();

    if (gmod / g_max > rr)  // g == gmod?
    {
      double r1 = random.Uniform();
      double r2 = random.Uniform();

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
      max_speed_squared = std::max({max_speed_squared, speed_squared(v1), speed_squared(v2)});
      g_max = 2 * std::sqrt(max_speed_squared);
      frequency_t = factor / g_max;
    }

  }
}

}
