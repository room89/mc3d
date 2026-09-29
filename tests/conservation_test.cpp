#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "cell.h"
#include "mirror_boundary.h"
#include "solver_setup.h"
#include "utils/utils.hpp"

namespace {

struct Totals {
  std::size_t count = 0;
  mc3d::Point momentum{0, 0, 0};
  double energy = 0;
};

Totals Measure(mc3d::Cell& cell) {
  const auto count = cell.GetParticleCount();
  return {count, cell.GetVelocity() * count, cell.GetEnergy() * count};
}

void ExpectConserved(const Totals& actual, const Totals& initial) {
  EXPECT_EQ(actual.count, initial.count);
  const double tolerance = 1e-10 * std::max(1.0, initial.energy);
  EXPECT_NEAR(actual.momentum.x, initial.momentum.x, tolerance);
  EXPECT_NEAR(actual.momentum.y, initial.momentum.y, tolerance);
  EXPECT_NEAR(actual.momentum.z, initial.momentum.z, tolerance);
  EXPECT_NEAR(actual.energy, initial.energy, tolerance);
}

}  // namespace

TEST(CellTimeStepTest, ThermalSpeedLimitsEveryAxisEvenWithZeroMeanVelocity) {
  for (int axis = 0; axis < 3; ++axis) {
    SCOPED_TRACE(axis);
    mc3d::Cell cell;
    cell.SetApex({0, 0, 0});
    const std::array<mc3d::Point, 3> sizes{
        {{1, 10, 10}, {10, 1, 10}, {10, 10, 1}}};
    cell.SetSize(sizes[axis]);
    // Six velocities give zero mean and T = 1 exactly.
    const std::array<mc3d::Point, 6> velocities{
        {{3, 0, 0}, {-3, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}}};
    for (const auto& velocity : velocities) {
      mc3d::Particle particle({0.5, 0.5, 0.5}, velocity);
      ASSERT_TRUE(cell.TryAcceptParticle(particle));
    }
    EXPECT_NEAR(cell.CalculateDt(), 1.0 / std::sqrt(2.0), 1e-14);
  }
}

TEST(CellTimeStepTest, BulkVelocityLimitsEachAxisIncludingNegativeVelocities) {
  for (const auto velocity :
       {mc3d::Point(-4, 0, 0), mc3d::Point(0, -6, 0), mc3d::Point(0, 0, -10)}) {
    mc3d::Cell cell;
    cell.SetApex({0, 0, 0});
    cell.SetSize(2, 3, 5);
    mc3d::Particle particle({1, 1, 1}, velocity);
    ASSERT_TRUE(cell.TryAcceptParticle(particle));
    EXPECT_DOUBLE_EQ(cell.CalculateDt(), 0.5);
  }
}

TEST(CellCollisionTest,
     RepeatedCollisionsConserveParticleCountMomentumAndEnergy) {
  utils::SeedRandom(20260929u);
  mc3d::Cell cell;
  cell.SetApex({0, 0, 0});
  cell.SetSize(1, 1, 1);
  cell.SetKn(0.05);
  cell.SetDt(0.02);
  for (int i = 0; i < 24; ++i) {
    mc3d::Particle particle(
        {0.5, 0.5, 0.5},
        {double(i % 5) - 1, double(i % 7) - 2, double(i % 3) - 0.5});
    ASSERT_TRUE(cell.TryAcceptParticle(particle));
  }
  const auto initial = Measure(cell);
  for (int step = 0; step < 100; ++step) {
    SCOPED_TRACE(step);
    cell.Collisions();
    ExpectConserved(Measure(cell), initial);
  }
}

TEST(PeriodicBoundaryTest,
     WrapsAllFacesAndMultipleLengthsWithoutChangingVelocity) {
  mc3d::SimulationConfig cfg;
  cfg.Lx = 2;
  cfg.Ly = 3;
  cfg.Lz = 5;
  cfg.apex_x = 4;
  cfg.apex_y = -7;
  cfg.apex_z = 2;
  const mc3d::Point velocity(1.25, -2.5, 3.75);
  const std::array<mc3d::Point, 8> positions{{{3.75, -6, 3},
                                              {6.25, -6, 3},
                                              {5, -7.25, 3},
                                              {5, -3.75, 3},
                                              {5, -6, 1.75},
                                              {5, -6, 7.25},
                                              {10.25, -13.25, 17.25},
                                              {6, -4, 7}}};
  const std::array<mc3d::Point, 8> expected{{{5.75, -6, 3},
                                             {4.25, -6, 3},
                                             {5, -4.25, 3},
                                             {5, -6.75, 3},
                                             {5, -6, 6.75},
                                             {5, -6, 2.25},
                                             {4.25, -4.25, 2.25},
                                             {4, -7, 2}}};
  std::vector<mc3d::Particle> particles;
  for (const auto& position : positions)
    particles.emplace_back(position, velocity);
  for (int face = 0; face < 6; ++face) {
    auto boundary =
        mc3d::MakeBoundary(mc3d::MakeBoundaryDescriptor(
                               static_cast<mc3d::BoundaryFace>(face), cfg),
                           mc3d::BoundaryType::Periodic, cfg);
    boundary->BoundaryCondition(particles, 1);
  }
  ASSERT_EQ(particles.size(), expected.size());
  for (std::size_t i = 0; i < particles.size(); ++i) {
    SCOPED_TRACE(i);
    EXPECT_DOUBLE_EQ(particles[i].position.x, expected[i].x);
    EXPECT_DOUBLE_EQ(particles[i].position.y, expected[i].y);
    EXPECT_DOUBLE_EQ(particles[i].position.z, expected[i].z);
    EXPECT_DOUBLE_EQ(particles[i].velocity.x, velocity.x);
    EXPECT_DOUBLE_EQ(particles[i].velocity.y, velocity.y);
    EXPECT_DOUBLE_EQ(particles[i].velocity.z, velocity.z);
  }
}

TEST(MirrorBoundaryTest, ExactUpperFaceHitRemainsInsideHalfOpenCell) {
  mc3d::Cell cell;
  cell.SetApex({-1, -1, -1});
  cell.SetSize(2, 2, 2);
  std::vector<mc3d::Particle> particles;
  particles.emplace_back(mc3d::Point(1, 0, 0), mc3d::Point(2, 3, 4));
  mc3d::MirrorBoundary boundary({1, 0, 0}, {1, 0, 0});
  boundary.BoundaryCondition(particles, 0.5);
  ASSERT_TRUE(cell.TryAcceptParticle(particles.front()));
  EXPECT_DOUBLE_EQ(cell.GetVelocity().x, -2);
  EXPECT_DOUBLE_EQ(cell.GetVelocity().y, 3);
  EXPECT_DOUBLE_EQ(cell.GetVelocity().z, 4);
  EXPECT_DOUBLE_EQ(cell.GetEnergy(), 14.5);
}
