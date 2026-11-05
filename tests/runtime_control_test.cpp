#include "runtime_control.h"
#include "config.h"
#include "utils/threadpool.hpp"
#include <gtest/gtest.h>
#include <atomic>
#include <filesystem>
#include <fstream>
#include <future>
#include <condition_variable>

namespace {
using namespace std::chrono_literals;
class RuntimeControlTest : public ::testing::Test {
 protected:
  std::filesystem::path directory = std::filesystem::temp_directory_path() /
      ("mc3d-control-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
  std::filesystem::path path = directory / "load.control";
  void SetUp() override { std::filesystem::create_directory(directory); }
  void TearDown() override { std::filesystem::remove_all(directory); }
  void Write(const std::string& mode) {
    auto temporary = directory / "new";
    { std::ofstream out(temporary); out << mode; }
    std::filesystem::rename(temporary, path);
  }
};

TEST_F(RuntimeControlTest, ModesAndInvalidOrMissingRequests) {
  mc3d::RuntimeControl control(path.string(), 5);
  control.Refresh();
  EXPECT_EQ(control.Threads(), 5);
  Write("balanced\n"); control.Refresh();
  EXPECT_EQ(control.Threads(), 3);
  EXPECT_DOUBLE_EQ(control.DutyCycle(), 0.6);
  Write("quiet"); control.Refresh();
  EXPECT_EQ(control.Threads(), 1);
  EXPECT_DOUBLE_EQ(control.DutyCycle(), 0.25);
  for (const auto& invalid : {"garbage", "quiet extra", "", "FULL"}) {
    Write(invalid); control.Refresh(); EXPECT_EQ(control.Mode(), "quiet");
  }
  Write("full" + std::string(100, ' ') + "garbage"); control.Refresh();
  EXPECT_EQ(control.Mode(), "quiet");
  std::filesystem::remove(path); control.Refresh();
  EXPECT_EQ(control.Mode(), "quiet");
  Write("full"); control.Refresh();
  EXPECT_EQ(control.Threads(), 5);
  EXPECT_DOUBLE_EQ(control.DutyCycle(), 1);
}

TEST_F(RuntimeControlTest, PauseAndThrottleCanBeReleasedDuringWait) {
  mc3d::RuntimeControl control(path.string(), 4);
  for (const auto& mode : {"pause", "quiet"}) {
    Write(mode);
    auto waiting = std::async(std::launch::async, [&] { return control.WaitBeforeStep(1s); });
    EXPECT_EQ(waiting.wait_for(250ms), std::future_status::timeout);
    Write("full");
    EXPECT_EQ(waiting.wait_for(2s), std::future_status::ready);
    EXPECT_EQ(waiting.get(), 4);
  }
}

TEST_F(RuntimeControlTest, QuietAddsRestBetweenSteps) {
  Write("quiet");
  mc3d::RuntimeControl control(path.string(), 1);
  const auto started = std::chrono::steady_clock::now();
  EXPECT_EQ(control.WaitBeforeStep(30ms), 1);
  EXPECT_GE(std::chrono::steady_clock::now() - started, 90ms);
}

TEST(ThreadPoolControlTest, LimitsConcurrencyThenRestoresCapacity) {
  utils::ThreadPool pool(4, "test");
  EXPECT_THROW(pool.SetActiveThreadCount(0), std::invalid_argument);
  EXPECT_THROW(pool.SetActiveThreadCount(5), std::invalid_argument);
  for (const auto workers : {1u, 4u, 2u}) {
    pool.SetActiveThreadCount(workers);
    EXPECT_EQ(pool.GetCapacity(), 4);
    std::mutex mutex;
    std::condition_variable changed;
    unsigned entered = 0;
    bool release = false;
    std::vector<std::future<void>> futures;
    for (int i = 0; i < 8; ++i) futures.push_back(pool.Execute([&] {
      std::unique_lock lock(mutex);
      ++entered;
      changed.notify_all();
      changed.wait(lock, [&] { return release; });
    }));
    {
      std::unique_lock lock(mutex);
      EXPECT_TRUE(changed.wait_for(lock, 2s, [&] { return entered >= workers; }));
      EXPECT_FALSE(changed.wait_for(lock, 50ms, [&] { return entered > workers; }));
      release = true;
    }
    changed.notify_all();
    for (auto& future : futures) future.get();
    EXPECT_EQ(entered, 8);
  }
}

TEST(ThreadPoolControlTest, StopDrainsQueueWithInactiveWorkers) {
  utils::ThreadPool pool(4, "test");
  pool.SetActiveThreadCount(1);
  std::atomic<int> count{0};
  for (int i = 0; i < 100; ++i) pool.Execute([&] { ++count; });
  pool.Stop();
  EXPECT_EQ(count, 100);
}

TEST(SimulationConfigTest, RuntimeControlPath) {
  char program[] = "solver", option[] = "--runtime-control", value[] = "load.control";
  char* args[] = {program, option, value};
  EXPECT_EQ(mc3d::LoadSimulationConfig(3, args).runtime_control, "load.control");
}
}  // namespace
