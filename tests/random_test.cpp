#include <gtest/gtest.h>
#include "utils/utils.hpp"
#include "utils/threadpool.hpp"

TEST(RandomSeed, RepeatedSeedRestartsSequence) {
  utils::SeedRandom(123);
  const auto first = utils::RandomEngine()();
  utils::RandomEngine()();
  utils::SeedRandom(123);
  EXPECT_EQ(first, utils::RandomEngine()());
  utils::SeedRandom(456);
  EXPECT_NE(first, utils::RandomEngine()());
}

TEST(RandomSeed, ExistingWorkerRestartsSequence) {
  utils::ThreadPool pool(1, "seed-test");
  uint32_t first = 0, second = 0;
  utils::SeedRandom(123);
  pool.Execute([&] { first = utils::RandomEngine()(); }).get();
  utils::SeedRandom(123);
  pool.Execute([&] { second = utils::RandomEngine()(); }).get();
  EXPECT_EQ(first, second);
}
