#include "particle_sampler.h"
#include "collision_model.h"
#include <gtest/gtest.h>
#include "utils/utils.hpp"
#include "utils/threadpool.hpp"
TEST(RandomStreamTest, RepeatedSeedRestartsSequence) {
  utils::RandomStream random(123);
  const auto first = random.Engine()();
  random.Engine()();
  random.Seed(123);
  EXPECT_EQ(first, random.Engine()());
  random.Seed(456);
  EXPECT_NE(first, random.Engine()());
}
TEST(RandomStreamTest, ExistingWorkerUsesExplicitStream) {
  utils::ThreadPool pool(1, "seed-test");
  utils::RandomStream random(123);
  uint32_t first = 0, second = 0;
  pool.Execute([&] { first = random.Engine()(); }).get();
  random.Seed(123);
  pool.Execute([&] { second = random.Engine()(); }).get();
  EXPECT_EQ(first, second);
}
TEST(RandomStreamsTest, CellIdsAreIndependentOfAccessOrder) {
  utils::RandomStreams a(42), b(42);
  auto a0 = a.At(0), a7 = a.At(7);
  auto b7 = b.At(7), b0 = b.At(0);
  EXPECT_EQ(a0->Engine()(), b0->Engine()());
  EXPECT_EQ(a7->Engine()(), b7->Engine()());
  EXPECT_EQ(a.At(0), a0);
  EXPECT_NE(a0->Engine()(), a7->Engine()());
}

TEST(PhysicsKernelsTest, SamplingAndCollisionsWorkWithoutCells) {
  utils::RandomStream random(57), repeated(57);
  std::vector<mc3d::Particle> particles, copy;
  mc3d::ParticleSampler::Initial(particles, {0,0,0}, {1,1,1}, 64, 1, {2,0,0}, random);
  mc3d::ParticleSampler::Initial(copy, {0,0,0}, {1,1,1}, 64, 1, {2,0,0}, repeated);
  ASSERT_EQ(particles.size(), 64u);
  mc3d::Point momentum(0,0,0);
  double energy = 0;
  for (size_t i = 0; i < particles.size(); ++i) {
    EXPECT_DOUBLE_EQ(particles[i].velocity.x, copy[i].velocity.x);
    EXPECT_DOUBLE_EQ(particles[i].position.z, copy[i].position.z);
    momentum += particles[i].velocity;
    energy += particles[i].velocity * particles[i].velocity;
  }
  EXPECT_NEAR(momentum.x / 64, 2, 1e-12);
  EXPECT_NEAR((energy / 64 - 4) / 3, 1, 1e-12);
  mc3d::CollisionModel::Apply(particles, momentum / 64, 1, 0.05, 0.2, random);
  mc3d::Point after(0,0,0);
  double energy_after = 0;
  for (const auto& particle : particles) {
    after += particle.velocity;
    energy_after += particle.velocity * particle.velocity;
  }
  EXPECT_EQ(particles.size(), 64u);
  EXPECT_NEAR((after - momentum).Mod(), 0, 1e-10);
  EXPECT_NEAR(energy_after, energy, 1e-10);
}
