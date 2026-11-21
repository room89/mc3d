#include "collision_scattering.h"

#include <gtest/gtest.h>

#include <array>

#include "collision_model.h"
using namespace mc3d;
namespace {
struct Moments {
  double a = 0, b = 0, c = 0;
  int n = 0;
  void Add(double x) {
    a += x;
    b += x * x;
    c += x * x * x * x;
    ++n;
  }
  void Check(double extra = 0) {
    ASSERT_GT(n, 0);
    EXPECT_NEAR(a / n, 0, 6 * std::sqrt(1. / (3 * n)) + extra);
    EXPECT_NEAR(b / n, 1. / 3, 6 * std::sqrt(4. / (45 * n)) + extra);
    EXPECT_NEAR(c / n, 1. / 5, 6 * std::sqrt(16. / (225 * n)) + extra);
  }
};
const std::array<Point, 5> directions{{{1, 0, 0},
                                       {0, 1, 0},
                                       {0, 0, 1},
                                       Point(1, 1, 1).Normalize(),
                                       Point(1e-8, 1, 1e-8).Normalize()}};
}  // namespace
TEST(HsScatteringTest, EndpointsAndNorm) {
  for (double g : {0., 1., 10.})
    for (double u : {0., 1e-15, .5, 1. - 1e-15, 1.})
      for (double v : {0., .25, .5, 1.}) {
        auto p = detail::IsotropicRelativeVelocity(g, u, v);
        EXPECT_NEAR(p.Mod(), g, 1e-12 * std::max(1., g));
      }
}
TEST(HsScatteringTest, IndependentAngularMoments) {
  for (unsigned seed : {123, 456, 789}) {
    utils::RandomStream rng(seed);
    std::array<Moments, 5> m;
    double xy = 0, xz = 0, yz = 0;
    for (int i = 0; i < 100000; ++i) {
      double u = rng.Uniform(), v = rng.Uniform();
      auto p = detail::IsotropicRelativeVelocity(1, u, v);
      for (int j = 0; j < 5; ++j) m[j].Add(p * directions[j]);
      xy += p.x * p.y;
      xz += p.x * p.z;
      yz += p.y * p.z;
    }
    for (auto& x : m) x.Check();
    for (double x : {xy, xz, yz})
      EXPECT_NEAR(x / 100000, 0, 6 * std::sqrt(1. / 1500000));
  }
}
TEST(HsScatteringTest, PublicPairConservesAndScatters) {
  for (auto d : directions)
    for (auto drift : {Point(0, 0, 0), Point(3, -2, 1)}) {
      utils::RandomStream rng(123), repeated(123);
      Moments m;
      for (int i = 0; i < 100000; ++i) {
        Particle a({1, 2, 3}, drift - d), b({4, 5, 6}, drift + d);
        auto aa = a, bb = b;
        double bound = 2, freq = 1, other = 2, otherfreq = 1;
        EXPECT_TRUE(Collision(a, b, bound, freq, 1, rng));
        Collision(aa, bb, other, otherfreq, 1, repeated);
        EXPECT_DOUBLE_EQ(a.velocity.x, aa.velocity.x);
        EXPECT_DOUBLE_EQ(b.velocity.z, bb.velocity.z);
        EXPECT_NEAR((a.velocity + b.velocity - 2 * drift).Mod(), 0, 1e-12);
        EXPECT_NEAR(a.velocity * a.velocity + b.velocity * b.velocity,
                    2 * (drift * drift + 1),
                    1e-12 * std::max(1., 2 * (drift * drift + 1)));
        EXPECT_DOUBLE_EQ(a.position.x, 1);
        EXPECT_DOUBLE_EQ(b.position.z, 6);
        m.Add(((b.velocity - a.velocity) * .5) * d);
      }
      m.Check();
    }
}
TEST(HsScatteringTest, SolverPathAngularRegression) {
  for (int axis = 0; axis < 3; ++axis) {
    utils::RandomStream rng(123 + axis);
    Moments m;
    double max_momentum_error = 0, max_energy_error = 0;
    std::vector<Particle> p(6);
    auto d = directions[axis];
    for (int trial = 0; trial < 2000000; ++trial) {
      for (int i = 0; i < 6; ++i) p[i].velocity = d * (i < 3 ? 1. : -1.);
      CollisionModel::Apply(p, {0, 0, 0}, 1, .1, 1e-4, rng);
      int changed = 0, last = 0;
      Point momentum(0, 0, 0);
      double energy = 0;
      for (int i = 0; i < 6; ++i) {
        auto v = p[i].velocity;
        momentum += v;
        energy += v * v;
        if ((v - d * (i < 3 ? 1. : -1.)).Mod() > 1e-12) {
          ++changed;
          last = i;
        }
      }
      if (changed) {
        max_momentum_error = std::max(max_momentum_error, momentum.Mod());
        max_energy_error = std::max(max_energy_error, std::abs(energy - 6));
        EXPECT_LT(momentum.Mod(), 1e-10);
        EXPECT_NEAR(energy, 6, 6e-10);
      }
      if (changed == 2) m.Add((p[last].velocity * d) / p[last].velocity.Mod());
    }
    ASSERT_GE(m.n, 2000);
    // Two changed particles can include rare repeated impacts of the same pair.
    // This integration regression allows 0.005 contamination; the independent
    // sampler and forced pair tests above use the exact isotropic variances.
    m.Check(.005);
    RecordProperty("momentum_error_" + std::to_string(axis),
                   max_momentum_error);
    RecordProperty("energy_error_" + std::to_string(axis), max_energy_error);
    RecordProperty("samples_" + std::to_string(axis), m.n);
    RecordProperty("second_" + std::to_string(axis), m.b / m.n);
    RecordProperty("fourth_" + std::to_string(axis), m.c / m.n);
  }
}
