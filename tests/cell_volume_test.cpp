#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "cell.h"
#include "geometry.h"
#include "solver_setup.h"
using namespace mc3d;
TEST(CellVolumeTest, CubeAnalyticAndThinGas) {
  for (double shift : {-.3, 0., .17, .6, 1.}) {
    Geometry body;
    body.CreateCube(shift, .7, .8, .6);
    for (double x : {-.5, 0., .5, 1.})
      for (double y : {-.5, 0., .5})
        for (double z : {-.5, 0., .5}) {
          Point a(x, y, z), size(.4, .3, .2);
          auto r = body.GasVolume(a, size);
          auto overlap = [](double a, double b, double c, double d) {
            return std::max(0., std::min(b, d) - std::max(a, c));
          };
          double solid = overlap(x, x + .4, shift, shift + .8) *
                         overlap(y, y + .3, -.3, .3) *
                         overlap(z, z + .2, -.35, .35);
          EXPECT_NEAR(r.solid, solid, 1e-12 * r.full);
          EXPECT_NEAR(r.gas, r.full - solid, 1e-12 * r.full);
        }
  }
  for (double gap : {1e-4, 1e-6, 1e-8, 1e-10}) {
    Geometry body;
    body.CreateCube(gap, 4, 2, 4);
    auto r = body.GasVolume({0, -.5, -.5}, {1, 1, 1});
    EXPECT_GT(r.gas, 0);
    EXPECT_NEAR(r.gas, gap, 1e-3 * gap);
  }
  Geometry full;
  full.CreateCube(-2, 4, 4, 4);
  EXPECT_EQ(full.GasVolume({0, 0, 0}, {1, 1, 1}).gas, 0);
}
TEST(CellVolumeTest, WedgeSlicesAndConservation) {
  const double angle = .3, h = std::sin(angle);
  Geometry body;
  body.CreateWedge(0, .8, 1, angle);
  for (int n : {4, 9, 20}) {
    double total = 0;
    for (int i = 0; i < n; ++i) {
      double a = double(i) / n, b = double(i + 1) / n;
      auto r = body.GasVolume({a, -1, -1}, {b - a, 2, 2});
      // Integral of 2*x*sin(alpha) across the x slab, times prism width.
      EXPECT_NEAR(r.solid, .8 * h * (b * b - a * a), 1e-12 * r.full);
      total += r.solid;
    }
    EXPECT_NEAR(total, .8 * h, 1e-10 * .8 * h);
  }
  auto before = body.GasVolume({-.2, -.4, -.5}, {.7, .7, .9});
  body.Fragment(.2);
  auto after = body.GasVolume({-.2, -.4, -.5}, {.7, .7, .9});
  EXPECT_NEAR(before.gas, after.gas, 1e-12 * before.full);
  body.Move({.3, -.1, .5});
  body.Scale(2);
  auto moved = body.GasVolume({.2, -1, 0}, {1.4, 1.4, 1.8});
  EXPECT_NEAR(moved.gas, 8 * before.gas, 1e-12 * moved.full);
}
TEST(CellVolumeTest, SolidCellNoFallbackAndPopulationInvariant) {
  auto body = std::make_unique<Geometry>();
  body->CreateCube(-2, 4, 4, 4);
  Cell cell(123, 0);
  cell.SetApex({0, 0, 0});
  cell.SetSize(1, 1, 1);
  cell.SetParameters(0, 0, 1);
  cell.SetKn(.1);
  cell.Initialize(20, body);
  EXPECT_EQ(cell.GetVolume(), 0);
  EXPECT_EQ(cell.GetParticleCount(), 0);
  EXPECT_TRUE(std::isinf(cell.CalculateDt()));
  Particle p({.5, .5, .5}, {1, 0, 0});
  ASSERT_TRUE(cell.TryAcceptParticle(p));
  EXPECT_THROW(cell.Collisions(), std::runtime_error);
}
TEST(CellVolumeTest, GeometryPolicy) {
  SimulationConfig c;
  Geometry cube;
  cube.CreateCube(-.1, .2, .2, .2);
  c.geometry_type = "cube";
  EXPECT_NO_THROW(ValidateSimulationGeometry(c, cube));
  cube.Move({.5, 0, 0});
  EXPECT_THROW(ValidateSimulationGeometry(c, cube), std::invalid_argument);
  Geometry cylinder;
  cylinder.CreateCylinder(0, .1, .2, 12);
  EXPECT_THROW(cylinder.GasVolume({0, 0, 0}, {1, 1, 1}), std::invalid_argument);
  c.geometry_type = "cylinder";
  EXPECT_THROW(ValidateSimulationGeometry(c, cylinder), std::invalid_argument);
}
#include "collision_model.h"
#include "inner_boundary.h"
namespace mc3d {
class CellVolumeTestAccess {
 public:
  static auto& Particles(Cell& c) { return c.particles; }
  static void Prepared(Cell& c) {
    c.RefreshMoments();
    c.CollisionsWithPreparedMoments();
  }
};
}  // namespace mc3d
TEST(CellVolumeTest, NoRandomDrawsAndIndependentOfBoundaryMark) {
  Geometry g;
  g.CreateCube(-.2, .4, .4, .4);
  InnerBoundary boundary(g, {0, 0, 0}, 1);
  utils::RandomStream a(123), b(123);
  Point center;
  auto v = boundary.CalcCellVolume({-.5, -.5, -.5}, {1, 1, 1}, a, &center);
  EXPECT_NEAR(v, 1 - .064, 1e-12);
  EXPECT_EQ(a.Engine()(), b.Engine()());
  auto body = std::make_unique<Geometry>();
  body->CreateCube(-.2, .4, .4, .4);
  Cell cell(123, 0);
  cell.SetApex({-.5, -.5, -.5});
  cell.SetSize(1, 1, 1);
  cell.SetParameters(0, 0, 1);
  cell.Initialize(0, body);
  cell.SetBodyMark(false);
  EXPECT_NEAR(cell.CalculateVolume(), 1 - .064, 1e-12);
  cell.SetApex({2, 2, 2});
  EXPECT_EQ(cell.GetVolume(), 1);
}
TEST(CellVolumeTest, VolumeFrequencyAndBothCellPaths) {
  for (double fraction : {1., .5, .01})
    for (double w : {.5, .75})
      for (unsigned seed : {123, 456, 789}) {
        auto body = std::make_unique<Geometry>();
        body->CreateCube(fraction, 4, 2, 4);
        auto random = std::make_shared<utils::RandomStream>(seed);
        Cell cell(random);
        cell.SetApex({0, -.5, -.5});
        cell.SetSize(1, 1, 1);
        cell.SetParameters(0, 0, 1);
        cell.SetKn(.1);
        cell.SetCharacteristicLength(1);
        cell.Initialize(0, body);
        cell.SetReferenceParticleCount(2);
        cell.collision_parameters = {
            CollisionKind::VHS, 6.6335e-26, 4e-10, 300, w, 300};
        auto& particles = CellVolumeTestAccess::Particles(cell);
        auto reset = [&] {
          particles = {Particle({fraction * .25, 0, 0}, {-1, 0, 0}),
                       Particle({fraction * .75, 0, 0}, {1, 0, 0})};
        };
        reset();
        EXPECT_NEAR(cell.CalculateKn(), .1 * fraction, 1e-14);
        const double power = 2 - 2 * w;
        const double moment = std::pow(4., power / 2) *
                              std::tgamma((power + 3) / 2) / std::tgamma(1.5);
        const double rate = std::sqrt(8 / std::acos(-1.)) *
                            std::pow(2., power) / (.1 * 2 * fraction * moment);
        double sum = 0;
        for (int trial = 0; trial < 20000; ++trial) {
          reset();
          sum += CollisionModel::ApplyWithStats(
                     particles, {0, 0, 0}, 1, cell.CalculateKn(), .5 / rate,
                     *random, cell.collision_parameters)
                     .accepted_collisions;
        }
        EXPECT_NEAR(sum / 20000, .5, 6 * std::sqrt(.5 / 20000));
        for (bool prepared : {false, true}) {
          reset();
          auto expected = particles;
          auto rng = *random;
          cell.SetDt(10 / rate);
          auto stats = CollisionModel::ApplyWithStats(
              expected, {0, 0, 0}, 1, cell.CalculateKn(), 10 / rate, rng,
              cell.collision_parameters);
          EXPECT_GT(stats.accepted_collisions, 0u);
          if (prepared)
            CellVolumeTestAccess::Prepared(cell);
          else
            cell.Collisions();
          for (int i = 0; i < 2; ++i)
            EXPECT_EQ((particles[i].velocity - expected[i].velocity).Mod(), 0);
          EXPECT_EQ(random->Engine()(), rng.Engine()());
        }
      }
}
TEST(CellVolumeTest, InitialPopulationUsesGasFractionWithoutChangingWeight) {
  auto body = std::make_unique<Geometry>();
  body->CreateCube(.25, 4, 2, 4);
  double count = 0;
  for (unsigned seed = 0; seed < 2000; ++seed) {
    Cell cell(seed, 0);
    cell.SetApex({0, -.5, -.5});
    cell.SetSize(1, 1, 1);
    cell.SetParameters(0, 0, 1);
    cell.SetKn(.1);
    cell.Initialize(40, body);
    cell.CleanInnerParticles(*body);
    count += cell.GetParticleCount();
    if (cell.GetParticleCount())
      EXPECT_NEAR(cell.CalculateKn(), .1 * 40 * .25 / cell.GetParticleCount(),
                  1e-14);
  }
  EXPECT_NEAR(count / 2000, 10, 6 * std::sqrt(40 * .25 * .75 / 2000));
}
TEST(CellVolumeTest, ShiftedGridTotalVolume) {
  for (int n : {5, 8, 13})
    for (double shift : {-.03, 0., .047})
      for (bool wedge : {false, true}) {
        Geometry body;
        if (wedge)
          body.CreateWedge(-.31, .44, .6, .4);
        else
          body.CreateCube(-.17, .45, .45, .45);
        body.Move({shift, .02, -.01});
        long double total = 0;
        int solid = 0, cut = 0;
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
              auto r = body.GasVolume({-.5 + double(i) / n, -.5 + double(j) / n,
                                       -.5 + double(k) / n},
                                      {1. / n, 1. / n, 1. / n});
              total += r.solid;
              solid += r.region == CellRegion::Solid;
              cut += r.region == CellRegion::Cut;
              EXPECT_GE(r.gas, 0);
              EXPECT_LE(r.gas, r.full * (1 + 1e-12));
            }
        double expected =
            wedge ? .44 * .6 * .6 * std::sin(.4) : .45 * .45 * .45;
        EXPECT_NEAR(double(total), expected, 1e-10 * expected);
        EXPECT_GT(cut, 0);
      }
}
TEST(CellVolumeTest, WedgeIndependentPiecewiseIntegral) {
  Geometry body;
  const double slope = std::sin(.3);
  body.CreateWedge(0, .8, 1, .3);
  for (double x : {-.2, .15, .5, .85})
    for (double y : {-.4, -.15, .05, .2})
      for (double z : {-.5, -.1, .3}) {
        Point size(.4, .2, .3);
        auto result = body.GasVolume({x, y, z}, size);
        double a = std::max(0., x), b = std::min(1., x + .4);
        std::vector<double> knots{a, b};
        for (double h : {y, y + .2, -y, -y - .2})
          if (h / slope > a && h / slope < b) knots.push_back(h / slope);
        std::sort(knots.begin(), knots.end());
        auto width = [&](double t) {
          return std::max(
              0., std::min(y + .2, slope * t) - std::max(y, -slope * t));
        };
        double area = 0;
        for (size_t i = 1; i < knots.size(); ++i)
          area += (knots[i] - knots[i - 1]) *
                  (width(knots[i]) + width(knots[i - 1])) / 2;
        double depth = std::max(0., std::min(.4, z + .3) - std::max(-.4, z));
        EXPECT_NEAR(result.solid, area * depth, 1e-12 * result.full);
      }
}

TEST(CellVolumeTest, CentroidTangencyAndNormals) {
  Geometry body;
  body.CreateCube(.25, 4, 2, 4);
  auto v = body.GasVolume({0, -.5, -.5}, {1, 1, 1});
  ASSERT_TRUE(v.gas_center);
  EXPECT_NEAR(v.gas_center->x, .125, 1e-14);
  EXPECT_NEAR(v.gas_center->y, 0, 1e-14);
  EXPECT_NEAR(v.gas_center->z, 0, 1e-14);
  body.ReverseNormals();
  EXPECT_EQ(body.GasVolume({0, -.5, -.5}, {1, 1, 1}).gas, v.gas);
  EXPECT_EQ(body.GasVolume({-.75, -.5, -.5}, {1, 1, 1}).region,
            CellRegion::Gas);
  auto solid = body.GasVolume({.25, -.5, -.5}, {1, 1, 1});
  EXPECT_EQ(solid.region, CellRegion::Solid);
  EXPECT_FALSE(solid.gas_center);
}
