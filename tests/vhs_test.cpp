#include <gtest/gtest.h>

#include "collision_model.h"
using namespace mc3d;
namespace {
CollisionParameters Params(double w) {
  CollisionParameters p;
  p.kind = CollisionKind::VHS;
  p.mass = 6.6335e-26;
  p.diameter = 4e-10;
  p.reference_temperature = 300;
  p.temperature_unit = 300;
  p.omega = w;
  return p;
}
}  // namespace
TEST(VhsTest, PairPoisson) {
  for (double w : {.5, .75, 1.})
    for (unsigned seed : {123, 456, 789})
      for (int steps : {1, 8, 32}) {
        auto par = Params(w);
        utils::RandomStream rng(seed);
        // Maxwell moment E[g^p] for unit temperature, independently normalized.
        double p = 2 - 2 * w;
        double rate = std::sqrt(8 / std::acos(-1.)) * std::pow(2., p) /
                      (2 * .1 * std::pow(4., p / 2) * std::tgamma((p + 3) / 2) /
                       std::tgamma(1.5));
        double sum = 0, sq = 0, zero = 0;
        for (int trial = 0; trial < 50000; ++trial) {
          std::vector<Particle> v{{{0, 0, 0}, {-1, 0, 0}},
                                  {{0, 0, 0}, {1, 0, 0}}};
          unsigned count = 0;
          for (int j = 0; j < steps; ++j)
            count += CollisionModel::ApplyWithStats(v, {0, 0, 0}, 1, .1,
                                                    1 / rate / steps, rng, par)
                         .accepted_collisions;
          sum += count;
          sq += double(count) * count;
          zero += count == 0;
        }
        EXPECT_NEAR(sum / 50000, 1, 6 / std::sqrt(50000.));
        EXPECT_NEAR((sq - sum * sum / 50000) / 49999, 1,
                    6 * std::sqrt(1 / 50000. + 2 / 49999.));
        EXPECT_NEAR(zero / 50000, std::exp(-1.),
                    6 * std::sqrt(std::exp(-1.) * (1 - std::exp(-1.)) / 50000));
      }
}
TEST(VhsTest, LimitsAndValidation) {
  auto par = Params(.5);
  for (double g : {0., .1, 1., 100.}) {
    EXPECT_DOUBLE_EQ(par.Rate(g), g);
    EXPECT_NEAR(par.PhysicalRate(g),
                std::acos(-1.) * par.diameter * par.diameter * g, 1e-30);
  }
  EXPECT_DOUBLE_EQ(Params(.75).Rate(0), 0);
  EXPECT_GT(Params(1).Rate(0), 0);
  par.omega = 1.1;
  EXPECT_THROW(par.Validate(), std::invalid_argument);
  par = Params(.75);
  par.mass = NAN;
  EXPECT_THROW(par.Validate(), std::invalid_argument);
  for (auto kind : {CollisionKind::None, CollisionKind::VHS})
    for (double dt : {0., -1.}) {
      par = Params(.75);
      par.kind = kind;
      utils::RandomStream a(123), b(123);
      std::vector<Particle> v{{{0, 0, 0}, {-1, 0, 0}}, {{0, 0, 0}, {1, 0, 0}}};
      EXPECT_EQ(CollisionModel::ApplyWithStats(v, {0, 0, 0}, 1, .1, dt, a, par)
                    .accepted_collisions,
                0);
      EXPECT_EQ(a.Engine()(), b.Engine()());
    }
}
TEST(VhsTest, MaxwellFrequencyAndEquilibrium) {
  for (double w : {.5, .75, 1.})
    for (int n : {2, 5, 32, 128})
      for (double temp : {1., 2., 4., 2000 / 300.})
        for (unsigned seed : {123, 456, 789}) {
          auto par = Params(w);
          utils::RandomStream rng(seed);
          std::normal_distribution<double> normal(0, std::sqrt(temp));
          const double rate = (n - 1) * std::sqrt(8 / std::acos(-1.)) /
                              (2 * .1) * std::pow(temp, 1 - w);
          double sum = 0, sq = 0, moment = 0, moment2 = 0;
          for (int trial = 0; trial < 20000; ++trial) {
            std::vector<Particle> v(n);
            Point mean(0, 0, 0);
            double energy = 0;
            for (auto& x : v) {
              x.velocity = {normal(rng.Engine()), normal(rng.Engine()),
                            normal(rng.Engine())};
              mean += x.velocity;
              energy += x.velocity * x.velocity;
            }
            auto before = mean;
            mean /= n;
            double c = CollisionModel::ApplyWithStats(v, mean, 1, .1, 1 / rate,
                                                      rng, par)
                           .accepted_collisions;
            sum += c;
            sq += c * c;
            Point after(0, 0, 0);
            double e = 0, xx = 0;
            for (auto& x : v) {
              after += x.velocity;
              e += x.velocity * x.velocity;
              xx += x.velocity.x * x.velocity.x;
            }
            EXPECT_NEAR((after - before).Mod(), 0,
                        1e-10 * std::max(1., energy));
            EXPECT_NEAR(e, energy, 1e-10 * std::max(1., energy));
            xx /= n;
            moment += xx;
            moment2 += xx * xx;
          }
          EXPECT_NEAR(sum / 20000, 1,
                      6 * std::sqrt((sq - sum * sum / 20000) / 19999 / 20000));
          EXPECT_NEAR(moment / 20000, temp,
                      6 * std::sqrt((moment2 - moment * moment / 20000) /
                                    19999 / 20000));
        }
}
namespace {
// Independent direct Gillespie oracle: enumerate all unordered pairs.
double Beam(bool direct, unsigned seed, double w, int steps, bool reverse) {
  utils::RandomStream rng(seed);
  auto par = Params(w);
  const int n = 16;
  std::vector<Particle> v(n);
  for (int i = 0; i < n; ++i) v[i].velocity = {i < n / 2 ? -1. : 1., 0, 0};
  if (reverse) std::reverse(v.begin(), v.end());
  const double power = 2 - 2 * w;
  const double coefficient = std::sqrt(8 / std::acos(-1.)) /
                             (.1 * n * std::pow(4., power / 2) *
                              std::tgamma((power + 3) / 2) / std::tgamma(1.5));
  if (!direct) {
    for (int s = 0; s < steps; ++s)
      CollisionModel::ApplyWithStats(v, {0, 0, 0}, 1, .1, .3 / steps, rng, par);
  } else {
    double t = 0;
    while (true) {
      double rates[n][n] = {}, total = 0;
      for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
          total += rates[i][j] =
              coefficient *
              std::pow((v[i].velocity - v[j].velocity).Mod(), power);
      if (total == 0) break;
      t += -std::log(std::max(rng.Uniform(), 1e-300)) / total;
      if (t > .3) break;
      double choice = rng.Uniform() * total;
      int a = 0, b = 1;
      bool found = false;
      for (int i = 0; i < n && !found; ++i)
        for (int j = i + 1; j < n; ++j)
          if ((choice -= rates[i][j]) <= 0) {
            a = i;
            b = j;
            found = true;
            break;
          }
      const double g = (v[a].velocity - v[b].velocity).Mod();
      const Point mean = (v[a].velocity + v[b].velocity) * .5;
      double z = 2 * rng.Uniform() - 1,
             phi = 2 * std::acos(-1.) * rng.Uniform();
      Point relative(g * std::sqrt(1 - z * z) * std::cos(phi),
                     g * std::sqrt(1 - z * z) * std::sin(phi), g * z);
      v[a].velocity = mean + relative * .5;
      v[b].velocity = mean - relative * .5;
    }
  }
  double xx = 0, e = 0;
  Point momentum(0, 0, 0);
  for (auto& x : v) {
    xx += x.velocity.x * x.velocity.x;
    e += x.velocity * x.velocity;
    momentum += x.velocity;
  }
  EXPECT_NEAR(e, n, 1e-9);
  EXPECT_NEAR(momentum.Mod(), 0, 1e-9);
  return xx / n;
}
}  // namespace
TEST(VhsTest, BeamsAgainstDirectPairOracle) {
  for (double w : {.5, .75, 1.})
    for (unsigned seed : {123, 456, 789})
      for (int steps : {1, 8})
        for (bool reverse : {false, true}) {
          double sum = 0, sq = 0;
          for (int trial = 0; trial < 2000; ++trial) {
            double delta = Beam(false, seed + trial * 17, w, steps, reverse) -
                           Beam(true, seed + trial * 19, w, 1, false);
            sum += delta;
            sq += delta * delta;
          }
          EXPECT_NEAR(sum / 2000, 0,
                      6 * std::sqrt((sq - sum * sum / 2000) / 1999 / 2000));
        }
}
TEST(VhsTest, PhysicalGammaAndNone) {
  auto p = Params(.75);
  double g = std::sqrt(4 * 1.380649e-23 * 300 / p.mass);
  EXPECT_NEAR(
      p.PhysicalRate(g) / (std::acos(-1.) * p.diameter * p.diameter * g),
      1 / 0.9190625268488832, 1e-14);
  p.kind = CollisionKind::None;
  utils::RandomStream a(123), b(123);
  std::vector<Particle> v{{{0, 0, 0}, {-1, 0, 0}}, {{0, 0, 0}, {1, 0, 0}}};
  EXPECT_EQ(CollisionModel::ApplyWithStats(v, {0, 0, 0}, 1, .1, 1, a, p)
                .accepted_collisions,
            0);
  EXPECT_EQ(a.Engine()(), b.Engine()());
}
