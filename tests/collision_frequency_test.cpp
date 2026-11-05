#include <gtest/gtest.h>

#include <cmath>
#include <numbers>

#include "collision_model.h"
using namespace mc3d;
namespace {
struct Samples {
  double sum = 0, squared = 0;
  int zero = 0, n = 0;
  void Add(double x) {
    sum += x;
    squared += x * x;
    zero += x == 0;
    ++n;
  }
  double Mean() const { return sum / n; }
  double Variance() const { return (squared - sum * sum / n) / (n - 1); }
  void Poisson(double mu) {
    EXPECT_NEAR(Mean(), mu, 6 * std::sqrt(mu / n));
    EXPECT_NEAR(Variance(), mu, 6 * std::sqrt(mu / n + 2 * mu * mu / (n - 1)));
    double p = std::exp(-mu);
    EXPECT_NEAR(double(zero) / n, p, 6 * std::sqrt(p * (1 - p) / n));
  }
};
// HS: lambda_ref=L*Kn=1/(sqrt(2)*n_ref*sigma).
// Each distinct pair has rate w*sigma*g/Vgas = g/(sqrt(2)*L*Kn*N0*f).
double PairRate(double g, double L, double kn, double n0, double f) {
  return g / (std::sqrt(2.) * L * kn * n0 * f);
}
}  // namespace
TEST(CollisionFrequencyTest, PairCountsArePoisson) {
  for (unsigned seed : {123, 456, 789})
    for (double mu : {.1, 1., 3.})
      for (int steps : {1, 8, 32}) {
        utils::RandomStream rng(seed);
        Samples samples;
        for (int trial = 0; trial < 50000; ++trial) {
          std::vector<Particle> p{{{0, 0, 0}, {-1, 0, 0}},
                                  {{0, 0, 0}, {1, 0, 0}}};
          double count = 0;
          for (int step = 0; step < steps; ++step)
            count += CollisionModel::ApplyWithStats(
                         p, {0, 0, 0}, 1, .1,
                         mu / PairRate(2, 1, .1, 2, 1) / steps, rng)
                         .accepted_collisions;
          samples.Add(count);
        }
        samples.Poisson(mu);
      }
}
TEST(CollisionFrequencyTest, PhysicalScalesAndDrift) {
  for (double L : {.5, 2.})
    for (double kn : {.05, .2})
      for (double n0 : {2., 20.})
        for (double f : {.25, 1.}) {
          utils::RandomStream rng(456);
          Samples samples;
          for (int trial = 0; trial < 50000; ++trial) {
            std::vector<Particle> p{{{0, 0, 0}, {2, -2, 1}},
                                    {{0, 0, 0}, {4, -2, 1}}};
            samples.Add(CollisionModel::ApplyWithStats(
                            p, {3, -2, 1}, L, kn * n0 * f / 2,
                            1 / PairRate(2, L, kn, n0, f), rng)
                            .accepted_collisions);
          }
          samples.Poisson(1);
        }
}
TEST(CollisionFrequencyTest, MaxwellianAbsoluteEnsembleRate) {
  for (int n : {2, 3, 4, 5, 6, 32, 128})
    for (unsigned seed : {123, 456, 789}) {
      utils::RandomStream rng(seed);
      std::normal_distribution<double> normal(0, 1);
      Samples samples;
      // Independent Maxwellian velocities, without sample temperature
      // rescaling.
      const double rate = n * (n - 1) / 2. * 4 / std::sqrt(std::numbers::pi) /
                          (std::sqrt(2.) * .1 * n);
      for (int trial = 0; trial < 20000; ++trial) {
        std::vector<Particle> p(n);
        Point mean(0, 0, 0);
        for (auto& x : p) {
          x.velocity = {normal(rng.Engine()), normal(rng.Engine()),
                        normal(rng.Engine())};
          mean += x.velocity;
        }
        mean /= n;
        samples.Add(
            CollisionModel::ApplyWithStats(p, mean, 1, .1, 1 / rate, rng)
                .accepted_collisions);
      }
      EXPECT_NEAR(samples.Mean(), 1,
                  6 * std::sqrt(samples.Variance() / samples.n));
      RecordProperty("mean_" + std::to_string(n) + "_" + std::to_string(seed),
                     samples.Mean());
    }
}
TEST(CollisionFrequencyTest, NoEventsAndNoRandomDraws) {
  for (int n : {0, 1, 2, 5, 6})
    for (double dt : {-1., 0., 1.}) {
      utils::RandomStream rng(123), reference(123);
      std::vector<Particle> p(n);
      for (auto& x : p) x.velocity = {3, 2, 1};
      EXPECT_EQ(CollisionModel::ApplyWithStats(p, {3, 2, 1}, 1, .1, dt, rng)
                    .accepted_collisions,
                0u);
      EXPECT_EQ(rng.Engine()(), reference.Engine()());
    }
}
TEST(CollisionFrequencyTest, StatisticsPreserveTrajectoryAndInvariants) {
  for (int n : {2, 3, 4, 5, 6, 32, 128}) {
    utils::RandomStream rng(123), reference(123);
    std::vector<Particle> p(n);
    Point before(0, 0, 0);
    double energy = 0;
    for (int i = 0; i < n; ++i) {
      p[i].velocity = {double(i % 3) - 1, double(i % 5) - 2, double(i % 7) - 3};
      before += p[i].velocity;
      energy += p[i].velocity * p[i].velocity;
    }
    auto copy = p;
    std::uint64_t events = 0;
    for (int step = 0; step < 100; ++step) {
      events +=
          CollisionModel::ApplyWithStats(p, before / double(n), 1, .1, .02, rng)
              .accepted_collisions;
      CollisionModel::Apply(copy, before / double(n), 1, .1, .02, reference);
      Point momentum(0, 0, 0);
      double after = 0;
      for (int i = 0; i < n; ++i) {
        momentum += p[i].velocity;
        after += p[i].velocity * p[i].velocity;
        EXPECT_DOUBLE_EQ((p[i].velocity - copy[i].velocity).Mod(), 0);
        EXPECT_DOUBLE_EQ(p[i].position.Mod(), 0);
      }
      EXPECT_NEAR((momentum - before).Mod(), 0, 1e-10 * std::max(1., energy));
      EXPECT_NEAR(after, energy, 1e-10 * std::max(1., energy));
    }
    EXPECT_GT(events, 0u);
    EXPECT_EQ(rng.Engine()(), reference.Engine()());
  }
}
