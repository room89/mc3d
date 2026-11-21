#include "statistical_moments.h"

#include <gtest/gtest.h>

#include <iomanip>
#include <random>
#include <sstream>
using namespace mc3d;
TEST(StatisticalMomentsTest, WeightedCentralMomentsAndGalileanInvariance) {
  for (double shift : {0., 1e6}) {
    WindowMoments a;
    a.duration = 3;
    for (int i = 0; i < 3; ++i)
      a.Add({shift + (i == 2 ? 3. : 0.), 2, 0}, i == 2 ? 2 : 1);
    auto r = a.Evaluate(.5, 2);
    EXPECT_NEAR(r.velocity.x, shift + 1.5, 1e-10);
    EXPECT_NEAR(r.temperature, .75, 1e-12);
    EXPECT_NEAR(r.pressure[0], .75, 1e-12);
    EXPECT_NEAR(r.heat.x, 0, 1e-12);
    EXPECT_EQ(r.pressure[4], 0);
  }
  WindowMoments b;
  b.duration = 1;
  b.Add({0, 0, 0}, 2);
  b.Add({3, 0, 0}, 1);
  auto r = b.Evaluate(1, 1);
  EXPECT_DOUBLE_EQ(r.temperature, 2. / 3);
  EXPECT_NEAR(r.heat.x, 3, 1e-12);
}
TEST(StatisticalMomentsTest,
     IndependentMaxwellSamplesRemoveSmallPopulationBias) {
  for (unsigned seed : {123, 456, 789})
    for (int n : {1, 2, 4, 20}) {
      std::mt19937 rng(seed);
      std::normal_distribution<double> normal(0, 1);
      WindowMoments a;
      double legacy = 0;
      for (int sample = 0; sample < 20000; ++sample) {
        a.duration += 1;
        Point sum{0, 0, 0};
        double squared = 0;
        for (int i = 0; i < n; ++i) {
          Point v(normal(rng), normal(rng), normal(rng));
          a.Add(v, 1);
          sum += v;
          squared += v * v;
        }
        auto mean = sum / double(n);
        legacy += (squared / n - mean * mean) / 3 / 20000;
      }
      const auto r = a.Evaluate(1, 1);
      const double tolerance = 5 * std::sqrt(2. / (3 * n * 20000)) + .005;
      const auto suffix = std::to_string(seed) + "_" + std::to_string(n);
      auto precise = [](double value) {
        std::ostringstream s;
        s << std::setprecision(17) << value;
        return s.str();
      };
      RecordProperty("pooled_" + suffix, precise(r.temperature));
      RecordProperty("legacy_" + suffix, precise(legacy));
      RecordProperty("tolerance_" + suffix, precise(tolerance));
      EXPECT_NEAR(r.temperature, 1, tolerance);
      EXPECT_NEAR(legacy, double(n - 1) / n, tolerance);
    }
}
TEST(StatisticalMomentsTest, EmptyExposure) {
  WindowMoments a;
  a.duration = 2;
  EXPECT_EQ(a.a0, 0);
  EXPECT_FALSE(a.initialized);
  a.Add({1, 2, 3}, .25);
  EXPECT_DOUBLE_EQ(double(a.a0 / a.duration), .125);
  EXPECT_DOUBLE_EQ(a.Evaluate(1, 1).temperature, 0);
}
TEST(StatisticalMomentsTest, TensorAndHeatAgainstDirectCenteredSum) {
  const std::array<Point, 4> velocities{
      {{1, 2, 3}, {-2, 1, 5}, {4, -1, 0}, {0, 3, -2}}};
  const double weights[] = {.1, .2, .3, .4};
  Point mean(0, 0, 0);
  WindowMoments a;
  a.duration = 1;
  for (int i = 0; i < 4; ++i) {
    a.Add(velocities[i], weights[i]);
    mean += velocities[i] * weights[i];
  }
  const auto r = a.Evaluate(2, 3);
  Point q(0, 0, 0);
  double temperature = 0;
  std::array<double, 9> pressure{};
  for (int i = 0; i < 4; ++i) {
    auto c = velocities[i] - mean;
    double v[] = {c.x, c.y, c.z};
    temperature += weights[i] * (c * c) / 3;
    q += c * (weights[i] * (c * c) / 3);
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        pressure[3 * j + k] += 2. / 3 * weights[i] * v[j] * v[k];
  }
  for (int i = 0; i < 9; ++i) EXPECT_NEAR(r.pressure[i], pressure[i], 1e-12);
  EXPECT_NEAR(r.temperature, temperature, 1e-12);
  EXPECT_NEAR((r.heat - q).Mod(), 0, 1e-12);
}

#include "mirror_boundary.h"
#include "particle_sampler.h"
TEST(StatisticalMomentsTest, InflowEnergyMeasuredBeforeSameStepReflection) {
  Geometry body;
  body.CreateCube(.05, 10, .1, 10);
  body.PrepareQueries();
  InnerBoundary boundary(body, {.5, 0, 0}, 1);
  boundary.SetWallTemperature(.3);
  utils::RandomStream random(123);
  std::vector<Particle> arrivals;
  ConservedQuantity incoming, outgoing;
  ParticleSampler::HyperFree({0, -.5, -.5}, {1, 1, 1}, 1000, {10, 0, 0}, 1,
                             {-1, 0, 0}, .1, boundary, arrivals, random,
                             &incoming);
  for (const auto& p : arrivals) outgoing.Add(p.velocity);
  const auto impulse = body.TakeStepImpulse();
  const double heat = body.TakeStepHeat();
  ASSERT_EQ(incoming.count, 1000);
  ASSERT_EQ(outgoing.count, 1000);
  ASSERT_GT(heat, 1);
  EXPECT_NEAR(outgoing.energy + heat, incoming.energy, 1e-10 * incoming.energy);
  EXPECT_NEAR((outgoing.momentum + impulse - incoming.momentum).Mod(), 0,
              1e-10 * incoming.momentum.Mod());
}
TEST(StatisticalMomentsTest, MirrorTransfersMomentumWithoutHeat) {
  MirrorBoundary wall({1, 0, 0}, {1, 0, 0});
  std::vector<Particle> particles{Particle({1.1, 0, 0}, {2, 1, 0})};
  wall.ProcessOutgoing(particles, .1);
  EXPECT_DOUBLE_EQ(wall.step_impulse.x, 4);
  EXPECT_DOUBLE_EQ(wall.step_heat, 0);
  EXPECT_DOUBLE_EQ(particles[0].velocity.x, -2);
}
