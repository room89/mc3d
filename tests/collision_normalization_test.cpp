#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include "cell.h"
#include "utils/utils.hpp"

namespace mc3d {
class CellCollisionTestAccess {
 public:
  static void CalculatePrepared(Cell& cell) { cell.CalculateWithPreparedMoments(); }
  static double StoredTemperature(const Cell& cell) { return cell.T; }
  static std::vector<Particle>& Particles(Cell& cell) { return cell.particles; }
  static void GasVolume(Cell& cell, double volume) { cell.volume_ = volume; }
};
}

TEST(CollisionNormalizationTest, ReservoirCountDoesNotFollowPopulationChanges) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  mc3d::Cell cell(random_stream);
  cell.SetApex({0,0,0}); cell.SetSize(2,1,1);
  cell.SetParameters(0,0,1); cell.SetKn(.1);
  std::unique_ptr<mc3d::Geometry> body;
  cell.Initialize(64,body);
  EXPECT_NEAR(cell.CalculateKn(), .1, 1e-14);
  auto& particles=mc3d::CellCollisionTestAccess::Particles(cell);
  particles.resize(32);
  EXPECT_NEAR(cell.CalculateKn(), .2, 1e-14);
  mc3d::CellCollisionTestAccess::GasVolume(cell,1);
  EXPECT_NEAR(cell.CalculateKn(), .1, 1e-14);
  particles.clear();
  EXPECT_TRUE(std::isinf(cell.CalculateKn()));
  cell.SetReferenceParticleCount(32.5);
  EXPECT_THROW(cell.SetReferenceParticleCount(0),std::invalid_argument);
  EXPECT_THROW(cell.SetReferenceParticleCount(std::numeric_limits<double>::quiet_NaN()),std::invalid_argument);
  EXPECT_THROW(cell.SetReferenceParticleCount(-1),std::invalid_argument);
}

namespace {
// Homogeneous counterstreaming beams: collisions relax <vx^2> from 1 toward 1/3.
// Only particle weight changes between populations; no geometry or transport.
double Relaxation(int count, double kn, double duration, int steps=8, bool shuffle=true) {
  double result=0;
  constexpr int ensembles=400;
  for (int seed=0;seed<ensembles;++seed) {
    auto random_stream = std::make_shared<utils::RandomStream>(12000+seed);
    auto& random = *random_stream;
    mc3d::Cell cell(random_stream);
    cell.SetApex({0,0,0}); cell.SetSize(1,1,1);
    cell.SetParameters(0,0,1); cell.SetKn(kn);
    cell.SetCharacteristicLength(1);
    std::unique_ptr<mc3d::Geometry> body;
    cell.Initialize(count,body);
    auto& particles=mc3d::CellCollisionTestAccess::Particles(cell);
    for(int i=0;i<count;++i) particles[i].velocity={i<count/2 ? 1. : -1.,0,0};
    if (shuffle) std::shuffle(particles.begin(),particles.end(),random.Engine());
    cell.SetDt(duration/steps);
    for(int step=0;step<steps;++step) cell.Collisions();
    double xx=0, energy=0;
    mc3d::Point momentum(0,0,0);
    for(const auto& p:particles) {
      xx+=p.velocity.x*p.velocity.x;
      energy+=p.velocity*p.velocity;
      momentum+=p.velocity;
    }
    EXPECT_NEAR(energy/count,1,1e-12);
    EXPECT_LT(momentum.Mod(),1e-10);
    EXPECT_EQ(particles.size(),static_cast<std::size_t>(count));
    result+=xx/count;
  }
  return result/ensembles;
}
}

TEST(CollisionNormalizationTest, RelaxationIsIndependentOfModelParticleWeight) {
  const double baseline=Relaxation(128,.1,.2);
  EXPECT_GT(baseline,.45); // Guard against all cases spuriously saturating at equilibrium.
  EXPECT_LT(baseline,.9); // Guard against collisions being disabled.
  for(int count:{32,512}) {
    const double value=Relaxation(count,.1,.2);
    SCOPED_TRACE(count);
    EXPECT_NEAR(value,baseline,.035);
    RecordProperty("vx2_"+std::to_string(count),std::to_string(value));
  }
  RecordProperty("vx2_128",std::to_string(baseline));
  EXPECT_NEAR(Relaxation(128,.2,.4),baseline,.015);
}

TEST(CollisionNormalizationTest, RelaxationDoesNotDependOnArrayOrderOrStepSubdivision) {
  const double reference=Relaxation(128,.1,.2);
  EXPECT_NEAR(Relaxation(128,.1,.2,8,false),reference,.025);
  EXPECT_NEAR(Relaxation(128,.1,.2,1),reference,.025);
  EXPECT_NEAR(Relaxation(128,.1,.2,32),reference,.025);
}

TEST(CollisionNormalizationTest, ColdTranslatingGasDoesNotCollide) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  mc3d::Cell cell(random_stream);
  cell.SetApex({0,0,0}); cell.SetSize(1,1,1);
  cell.SetReferenceParticleCount(32); cell.SetKn(.1); cell.SetDt(1);
  for (int i=0;i<32;++i) {
    mc3d::Particle p({.5,.5,.5},{100,-20,3});
    ASSERT_TRUE(cell.TryAcceptParticle(p));
  }
  cell.Collisions();
  for(const auto& p:mc3d::CellCollisionTestAccess::Particles(cell)) {
    EXPECT_DOUBLE_EQ(p.velocity.x,100);
    EXPECT_DOUBLE_EQ(p.velocity.y,-20);
    EXPECT_DOUBLE_EQ(p.velocity.z,3);
  }
}

TEST(CellMomentsTest, EmptySingleAndCounterstreamingParticles) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  mc3d::Cell cell(random_stream);
  cell.SetTemperature(42);
  EXPECT_EQ(cell.ComputeMoments().count, 0u);
  EXPECT_DOUBLE_EQ(cell.ComputeMoments().temperature, 0);
  auto& particles = mc3d::CellCollisionTestAccess::Particles(cell);
  particles.emplace_back(mc3d::Point(0,0,0), mc3d::Point(3,0,0));
  auto moments = cell.ComputeMoments();
  EXPECT_EQ(moments.count, 1u);
  EXPECT_DOUBLE_EQ(moments.velocity.x, 3);
  EXPECT_DOUBLE_EQ(moments.energy, 4.5);
  EXPECT_DOUBLE_EQ(moments.temperature, 0);
  particles.emplace_back(mc3d::Point(0,0,0), mc3d::Point(-3,0,0));
  moments = cell.ComputeMoments();
  EXPECT_EQ(moments.count, 2u);
  EXPECT_DOUBLE_EQ(moments.velocity.x, 0);
  EXPECT_DOUBLE_EQ(moments.energy, 4.5);
  EXPECT_DOUBLE_EQ(moments.temperature, 3);
  EXPECT_DOUBLE_EQ(mc3d::CellCollisionTestAccess::StoredTemperature(cell), 42);
  EXPECT_DOUBLE_EQ(particles.front().velocity.x, 3);
}

TEST(CellMomentsTest, ResultsStayFreshAfterParticleChanges) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  mc3d::Cell cell(random_stream);
  cell.SetTemperature(42);
  auto& particles = mc3d::CellCollisionTestAccess::Particles(cell);
  particles.emplace_back(mc3d::Point(0,0,0), mc3d::Point(2,4,6));
  const auto before = cell.ComputeMoments();
  EXPECT_DOUBLE_EQ(before.velocity.x, 2);
  EXPECT_DOUBLE_EQ(before.velocity.y, 4);
  EXPECT_DOUBLE_EQ(before.velocity.z, 6);
  EXPECT_DOUBLE_EQ(before.energy, 28);
  particles.front().velocity = {-2,-4,-6};
  const auto after = cell.ComputeMoments();
  EXPECT_DOUBLE_EQ(after.velocity.x, -2);
  EXPECT_DOUBLE_EQ(after.velocity.y, -4);
  EXPECT_DOUBLE_EQ(after.velocity.z, -6);
  EXPECT_DOUBLE_EQ(before.velocity.x, 2);
  EXPECT_DOUBLE_EQ(mc3d::CellCollisionTestAccess::StoredTemperature(cell), 42);
  particles.clear();
  const auto empty = cell.ComputeMoments();
  EXPECT_DOUBLE_EQ(empty.velocity.Mod(), 0);
  EXPECT_DOUBLE_EQ(empty.energy, 0);
}

namespace {
void ExpectSameParticles(mc3d::Cell& actual, mc3d::Cell& expected) {
  const auto& a = mc3d::CellCollisionTestAccess::Particles(actual);
  const auto& b = mc3d::CellCollisionTestAccess::Particles(expected);
  ASSERT_EQ(a.size(), b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    SCOPED_TRACE(i);
    EXPECT_EQ(a[i].position.x, b[i].position.x);
    EXPECT_EQ(a[i].position.y, b[i].position.y);
    EXPECT_EQ(a[i].position.z, b[i].position.z);
    EXPECT_EQ(a[i].velocity.x, b[i].velocity.x);
    EXPECT_EQ(a[i].velocity.y, b[i].velocity.y);
    EXPECT_EQ(a[i].velocity.z, b[i].velocity.z);
  }
}

void PrepareCell(mc3d::Cell& cell, const std::unique_ptr<mc3d::Geometry>& body,
                 int count) {
  cell.SetApex({-1.5, -.5, -.5});
  cell.SetSize(1, 1, 1);
  cell.SetParameters(0, 0, 1);
  cell.SetKn(.1);
  cell.SetCharacteristicLength(1);
  ASSERT_TRUE(cell.Initialize(0, body));
  cell.SetReferenceParticleCount(std::max(count, 1));
  auto& particles = mc3d::CellCollisionTestAccess::Particles(cell);
  for (int i = 0; i < count; ++i)
    particles.emplace_back(mc3d::Point(-1, 0, 0),
                           mc3d::Point(2, (i % 3 - 1) * .2, (i % 5 - 2) * .1));
}
}

TEST(CellPreparedMomentsTest, MatchesStandaloneStepAndRandomState) {
  for (bool with_body : {false, true}) {
    std::unique_ptr<mc3d::Geometry> body;
    if (with_body) {
      body = std::make_unique<mc3d::Geometry>();
      body->CreateCube(-.5, 1, 1, 1);
      body->PrepareQueries();
    }
    for (int count : {0, 1, 5, 6, 64}) {
      SCOPED_TRACE(with_body);
      SCOPED_TRACE(count);
      auto rng_a = std::make_shared<utils::RandomStream>(123);
      auto rng_b = std::make_shared<utils::RandomStream>(123);
      mc3d::Cell prepared(rng_a), standalone(rng_b);
      PrepareCell(prepared, body, count);
      PrepareCell(standalone, body, count);
      EXPECT_EQ(prepared.CalculateDt(), standalone.CalculateDt());
      prepared.SetDt(.5);
      standalone.SetDt(.5);
      mc3d::CellCollisionTestAccess::CalculatePrepared(prepared);
      standalone.Calculate();
      ExpectSameParticles(prepared, standalone);
      EXPECT_EQ(rng_a->Engine()(), rng_b->Engine()());
      if (with_body && count == 1) {
        const auto& particle = mc3d::CellCollisionTestAccess::Particles(prepared).front();
        EXPECT_LT(particle.position.x, -.5);
        EXPECT_LT(particle.velocity.x, 0);
      }
    }
  }
}

TEST(CellPreparedMomentsTest, StandaloneCallsRefreshAfterPopulationAndVelocityChanges) {
  for (bool full_step : {false, true}) {
    auto rng_a = std::make_shared<utils::RandomStream>(321);
    auto rng_b = std::make_shared<utils::RandomStream>(321);
    mc3d::Cell actual(rng_a), expected(rng_b);
    std::unique_ptr<mc3d::Geometry> body;
    PrepareCell(actual, body, 64);
    PrepareCell(expected, body, 64);
    actual.CalculateDt();
    expected.CalculateDt();
    for (auto* cell : {&actual, &expected}) {
      auto& particles = mc3d::CellCollisionTestAccess::Particles(*cell);
      particles.resize(24);
      for (auto& particle : particles) particle.velocity = {7, -3, 2};
    }
    // Only the reference explicitly refreshes; public entry points must do it themselves.
    expected.CalculateDt();
    actual.SetDt(.02);
    expected.SetDt(.02);
    if (full_step) {
      actual.Calculate();
      mc3d::CellCollisionTestAccess::CalculatePrepared(expected);
    } else {
      actual.Collisions();
      expected.Collisions();
    }
    ExpectSameParticles(actual, expected);
    EXPECT_EQ(mc3d::CellCollisionTestAccess::StoredTemperature(actual), 0);
    EXPECT_EQ(rng_a->Engine()(), rng_b->Engine()());
  }
}

TEST(CellPreparedMomentsTest, PublicCalculateKeepsVirtualCollisionDispatch) {
  class CustomCell : public mc3d::Cell {
   public:
    void Collisions() override { collided = true; }
    void MoveParticles() override { moved = true; }
    bool collided = false;
    bool moved = false;
  } cell;
  cell.Calculate();
  EXPECT_TRUE(cell.collided);
  EXPECT_TRUE(cell.moved);
}
