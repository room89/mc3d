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

static_assert(std::uniform_random_bit_generator<utils::Pcg32>);
static_assert(sizeof(utils::RandomStream) == 16);

TEST(Pcg32Test, MatchesPinnedUpstreamVectors) {
  // Generated with pcg32_srandom_r/pcg32_random_r from pcg-c-basic
  // bc39cd76ac3d541e618606bcc6e1e5ba5e5e6aa3. First case is the upstream demo.
  struct Case { uint32_t seed; uint64_t id; uint32_t values[6]; };
  const Case cases[] = {
      {42, 54, {0xa15c02b7,0x7b47f409,0xba1d3330,0x83d2f293,0xbfa4784b,0xcbed606e}},
      {0, 0, {0xe4c14788,0x379c6516,0x5c4ab3bb,0x601d23e0,0x1c382b8c,0xd1faab16}},
      {UINT32_MAX, utils::Pcg32::kMaxStreamId,
       {0xc0b80fc6,0xfaf9a041,0xedab3fb5,0x13305841,0x5f00751f,0xf549787d}}};
  for (const auto& test : cases) {
    utils::RandomStream random(test.seed, test.id);
    for (int round = 0; round < 2; ++round) {
      random.Seed(test.seed, test.id);
      for (auto expected : test.values) EXPECT_EQ(random.Engine()(), expected);
    }
  }
}

TEST(Pcg32Test, RejectsAliasingStreamIdsWithoutChangingState) {
  utils::RandomStream actual(42, 54), expected(42, 54);
  EXPECT_EQ(actual.Engine()(), expected.Engine()());
  EXPECT_THROW(actual.Seed(1, uint64_t{1} << 63), std::invalid_argument);
  EXPECT_THROW(actual.Seed(1, UINT64_MAX), std::invalid_argument);
  EXPECT_EQ(actual.Engine()(), expected.Engine()());
  EXPECT_THROW((utils::RandomStream(1, UINT64_MAX)), std::invalid_argument);
}

TEST(RandomStreamsTest, InvalidIndexDoesNotDiscardExistingStreams) {
  utils::RandomStreams streams(42);
  const auto first = streams.At(0);
  EXPECT_THROW(streams.At(std::numeric_limits<std::size_t>::max()), std::exception);
  EXPECT_EQ(streams.At(0), first);
}
