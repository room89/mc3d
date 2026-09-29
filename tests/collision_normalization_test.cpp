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
  static std::vector<Particle>& Particles(Cell& cell) { return cell.particles; }
  static void GasVolume(Cell& cell, double volume) { cell.volume_ = volume; }
};
}

TEST(CollisionNormalizationTest, ReservoirCountDoesNotFollowPopulationChanges) {
  mc3d::Cell cell;
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
    utils::RandomEngine().seed(12000+seed);
    mc3d::Cell cell;
    cell.SetApex({0,0,0}); cell.SetSize(1,1,1);
    cell.SetParameters(0,0,1); cell.SetKn(kn);
    cell.SetCharacteristicLength(1);
    std::unique_ptr<mc3d::Geometry> body;
    cell.Initialize(count,body);
    auto& particles=mc3d::CellCollisionTestAccess::Particles(cell);
    for(int i=0;i<count;++i) particles[i].velocity={i<count/2 ? 1. : -1.,0,0};
    if (shuffle) std::shuffle(particles.begin(),particles.end(),utils::RandomEngine());
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
  mc3d::Cell cell;
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
