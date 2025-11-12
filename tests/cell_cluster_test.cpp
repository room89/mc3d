#include "cell_cluster.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "boundary.h"
#include "geometry.h"
#include "particle.h"
#include "point.h"

namespace mc3d {

class CellClusterTestAccess {
 public:
  static double GetDt(const CellCluster& cluster) { return cluster.dt; }

  static const std::vector<Cell*>& GetLookup(const CellCluster& cluster) {
    return cluster.cell_lookup_;
  }

  static bool FindCellIndex(const CellCluster& cluster, const Point& position,
                            std::size_t& index) {
    return cluster.FindCellIndex(position, index);
  }

  static std::size_t BufferSize(const CellCluster& cluster) {
    return cluster.partile_buffer.size();
  }

  static double GetSnapshotInterval(const CellCluster& cluster) {
    return cluster.snapshot_interval_;
  }

  static double GetNextSnapshotTime(const CellCluster& cluster) {
    return cluster.next_snapshot_time_;
  }

  static double GetTime(const CellCluster& cluster) { return cluster.t; }

  static bool RunTimeStep(CellCluster& cluster) { return cluster.TimeStep(); }

  static void DistributeParticles(CellCluster& cluster,
                                  std::vector<Particle>& buffer) {
    cluster.DistributeParticles(buffer);
  }
};

}  // namespace mc3d

namespace {

constexpr double kKn = 0.05;
constexpr double kCu = 0.2;
constexpr double kS = 10.0;
constexpr double kAlpha = 0.0;
constexpr double kTemperature = 1.0;

class ScopedWorkingDirectory {
 public:
  explicit ScopedWorkingDirectory(const std::filesystem::path& target)
      : original_(std::filesystem::current_path()), target_(target) {
    std::filesystem::remove_all(target_);
    std::filesystem::create_directories(target_);
    std::filesystem::current_path(target_);
  }

  ScopedWorkingDirectory(const ScopedWorkingDirectory&) = delete;
  ScopedWorkingDirectory& operator=(const ScopedWorkingDirectory&) = delete;

  ~ScopedWorkingDirectory() {
    std::error_code ec;
    std::filesystem::current_path(original_, ec);
    std::filesystem::remove_all(target_, ec);
  }

 private:
  std::filesystem::path original_;
  std::filesystem::path target_;
};

}  // namespace

TEST(CellClusterInitializeTest, SynchronizesDtAcrossCells) {
  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(2.0, 1.0, 1.0);
  cluster.SetEndTime(0.5);
  ASSERT_TRUE(cluster.Initialize(2, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  const double dt = mc3d::CellClusterTestAccess::GetDt(cluster);
  const auto& lookup = mc3d::CellClusterTestAccess::GetLookup(cluster);

  ASSERT_EQ(lookup.size(), 2U);
  EXPECT_GT(dt, 0.0);
  for (mc3d::Cell* cell_ptr : lookup) {
    ASSERT_NE(cell_ptr, nullptr);
    EXPECT_DOUBLE_EQ(cell_ptr->PeekDt(), dt);
  }
}

namespace mc3d::test {

class CountingBoundary : public mc3d::Boundary {
 public:
  int call_count = 0;
  double last_dt = 0.0;

  int BoundaryCondition(std::vector<mc3d::Particle>& /*particles*/,
                        double dt) override {
    ++call_count;
    last_dt = dt;
    return 0;
  }
};

}  // namespace mc3d::test

TEST(CellClusterInitializeTest, PopulatesLookupAndFindsCells) {
  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  const auto& lookup = mc3d::CellClusterTestAccess::GetLookup(cluster);
  ASSERT_EQ(lookup.size(), 1U);
  EXPECT_NE(lookup.front(), nullptr);

  std::size_t index = std::numeric_limits<std::size_t>::max();
  const bool inside = mc3d::CellClusterTestAccess::FindCellIndex(
      cluster, mc3d::Point(0.25, 0.25, 0.25), index);
  EXPECT_TRUE(inside);
  EXPECT_EQ(index, 0U);

  const bool outside = mc3d::CellClusterTestAccess::FindCellIndex(
      cluster, mc3d::Point(-0.1, 0.0, 0.0), index);
  EXPECT_FALSE(outside);
}

TEST(CellClusterInitializeTest, CleansParticleBufferAndMarksBodyCells) {
  auto geometry = std::make_unique<mc3d::Geometry>();
  geometry->CreateCube(-0.5, 1.0, 1.0, 1.0);

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(-0.5, -0.5, -0.5));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.5);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 200.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  const auto& lookup = mc3d::CellClusterTestAccess::GetLookup(cluster);
  ASSERT_EQ(lookup.size(), 1U);
  mc3d::Cell* cell = lookup.front();
  ASSERT_NE(cell, nullptr);

  EXPECT_EQ(mc3d::CellClusterTestAccess::BufferSize(cluster), 0U);
  EXPECT_TRUE(cell->GetBodyMark());
}

TEST(CellClusterGridTest, FindCellIndexRejectsBoundaryPoints) {
  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.1);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  std::size_t index = std::numeric_limits<std::size_t>::max();
  EXPECT_FALSE(mc3d::CellClusterTestAccess::FindCellIndex(
      cluster, mc3d::Point(0.0, 0.5, 0.5), index));
  EXPECT_FALSE(mc3d::CellClusterTestAccess::FindCellIndex(
      cluster, mc3d::Point(1.0, 0.5, 0.5), index));
  const double epsilon = 1e-6;
  EXPECT_TRUE(mc3d::CellClusterTestAccess::FindCellIndex(
      cluster, mc3d::Point(epsilon, epsilon, epsilon), index));
  EXPECT_EQ(index, 0U);
}

TEST(CellClusterDistributionTest, DistributeParticlesAcceptsAndRejects) {
  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(2.0, 1.0, 1.0);
  cluster.SetEndTime(0.1);
  ASSERT_TRUE(cluster.Initialize(2, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  std::vector<mc3d::Particle> buffer;
  buffer.emplace_back(mc3d::Point(0.25, 0.25, 0.25),
                      mc3d::Point(0.0, 0.0, 0.0));
  buffer.emplace_back(mc3d::Point(1.5, 0.25, 0.25), mc3d::Point(0.0, 0.0, 0.0));
  buffer.emplace_back(mc3d::Point(2.5, 0.25, 0.25), mc3d::Point(0.0, 0.0, 0.0));

  mc3d::CellClusterTestAccess::DistributeParticles(cluster, buffer);

  const auto& lookup = mc3d::CellClusterTestAccess::GetLookup(cluster);
  ASSERT_EQ(lookup.size(), 2U);
  ASSERT_NE(lookup[0], nullptr);
  ASSERT_NE(lookup[1], nullptr);

  EXPECT_EQ(lookup[0]->GetParticleCount(), 1U);
  EXPECT_EQ(lookup[1]->GetParticleCount(), 1U);
  ASSERT_EQ(buffer.size(), 1U);
  EXPECT_DOUBLE_EQ(buffer.front().position.x, 2.5);
}

TEST(CellClusterSnapshotTest, ZeroIntervalDisablesScheduling) {
  mc3d::CellCluster cluster;
  cluster.SetSnapshotInterval(0.0);

  EXPECT_DOUBLE_EQ(mc3d::CellClusterTestAccess::GetSnapshotInterval(cluster),
                   0.0);
  EXPECT_DOUBLE_EQ(mc3d::CellClusterTestAccess::GetNextSnapshotTime(cluster),
                   0.0);
}

TEST(CellClusterSnapshotTest, PositiveIntervalInitialisesNextTime) {
  mc3d::CellCluster cluster;
  cluster.SetSnapshotInterval(0.25);

  EXPECT_DOUBLE_EQ(mc3d::CellClusterTestAccess::GetSnapshotInterval(cluster),
                   0.25);
  EXPECT_DOUBLE_EQ(mc3d::CellClusterTestAccess::GetNextSnapshotTime(cluster),
                   0.25);
}

TEST(CellClusterSnapshotTest, TimeStepHonoursSnapshotSchedule) {
  const auto temp_dir =
      std::filesystem::temp_directory_path() / "mc3d_cluster_snapshot_test";
  ScopedWorkingDirectory guard(temp_dir);

  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.25);
  cluster.SetSnapshotInterval(0.1);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 25.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  bool finished = false;
  int steps = 0;
  while (!finished && steps < 500) {
    finished = mc3d::CellClusterTestAccess::RunTimeStep(cluster);
    ++steps;
  }
  EXPECT_TRUE(finished);
  EXPECT_GT(steps, 0);

  constexpr std::string_view kSnapshotPrefix = "data";
  constexpr std::string_view kSnapshotExtension = ".dat";
  std::vector<double> snapshot_times;
  for (const auto& entry : std::filesystem::directory_iterator(temp_dir)) {
    if (entry.is_regular_file()) {
      const auto filename = entry.path().filename().string();
      if (filename.rfind(kSnapshotPrefix, 0) == 0 &&
          filename.ends_with(kSnapshotExtension)) {
        const auto time_part = filename.substr(
            kSnapshotPrefix.size(), filename.size() - kSnapshotPrefix.size() -
                                        kSnapshotExtension.size());
        try {
          if (!time_part.empty()) {
            snapshot_times.push_back(std::stod(time_part));
          }
        } catch (const std::exception&) {
          // ignore parsing errors
        }
      }
    }
  }

  EXPECT_FALSE(snapshot_times.empty());
  std::sort(snapshot_times.begin(), snapshot_times.end());
  EXPECT_GE(snapshot_times.front(), 0.1);
  EXPECT_GE(snapshot_times.back(), 0.2);
  EXPECT_NEAR(mc3d::CellClusterTestAccess::GetNextSnapshotTime(cluster), 0.3,
              1e-6);
  EXPECT_NEAR(mc3d::CellClusterTestAccess::GetTime(cluster), 0.25, 1e-9);
}

TEST(CellClusterSnapshotWriters, WriteTextSnapshotProducesExpectedColumns) {
  const auto temp_dir =
      std::filesystem::temp_directory_path() / "mc3d_snapshot_text_test";
  ScopedWorkingDirectory guard(temp_dir);

  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.1);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));
  cluster.SetBinaryOutput(false);

  const auto snapshot_path = temp_dir / "snapshot.dat";
  ASSERT_TRUE(cluster.WriteFile(snapshot_path.string()));

  std::ifstream input(snapshot_path);
  ASSERT_TRUE(input.is_open());

  std::string header;
  ASSERT_TRUE(std::getline(input, header));
  EXPECT_EQ(header, "x;y;z;N;ro;T;vx;vy;vz;E");

  std::string row;
  ASSERT_TRUE(std::getline(input, row));
  std::vector<std::string> columns;
  std::stringstream ss(row);
  std::string cell;
  while (std::getline(ss, cell, ';')) {
    columns.push_back(cell);
  }
  ASSERT_EQ(columns.size(), 10U);

  const auto as_double = [](const std::string& value) {
    return std::stod(value);
  };

  EXPECT_NEAR(as_double(columns[0]), 0.5, 1e-9);
  EXPECT_NEAR(as_double(columns[1]), 0.5, 1e-9);
  EXPECT_NEAR(as_double(columns[2]), 0.5, 1e-9);
  EXPECT_EQ(columns[3], "0");
  EXPECT_NEAR(as_double(columns[4]), 0.0, 1e-12);
  EXPECT_NEAR(as_double(columns[5]), 0.0, 1e-12);
  EXPECT_NEAR(as_double(columns[6]), 0.0, 1e-12);
  EXPECT_NEAR(as_double(columns[7]), 0.0, 1e-12);
  EXPECT_NEAR(as_double(columns[8]), 0.0, 1e-12);
  EXPECT_NEAR(as_double(columns[9]), 0.0, 1e-12);

  EXPECT_FALSE(std::getline(input, row));
}

TEST(CellClusterSnapshotWriters, WriteBinarySnapshotProducesExpectedLayout) {
  const auto temp_dir =
      std::filesystem::temp_directory_path() / "mc3d_snapshot_binary_test";
  ScopedWorkingDirectory guard(temp_dir);

  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.1);
  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 0.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));
  cluster.SetBinaryOutput(true);

  const auto snapshot_path = temp_dir / "snapshot.bin";
  ASSERT_TRUE(cluster.WriteFile(snapshot_path.string()));

  std::ifstream input(snapshot_path, std::ios::binary);
  ASSERT_TRUE(input.is_open());

  char magic[4];
  input.read(magic, sizeof(magic));
  ASSERT_EQ(input.gcount(), static_cast<std::streamsize>(sizeof(magic)));
  EXPECT_EQ(std::string(magic, sizeof(magic)), "MC3D");

  std::uint32_t version = 0;
  input.read(reinterpret_cast<char*>(&version), sizeof(version));
  EXPECT_EQ(version, 1U);

  std::uint64_t record_count = 0;
  input.read(reinterpret_cast<char*>(&record_count), sizeof(record_count));
  EXPECT_EQ(record_count, 1U);

  double density_reference = -1.0;
  input.read(reinterpret_cast<char*>(&density_reference),
             sizeof(density_reference));
  EXPECT_DOUBLE_EQ(density_reference, 0.0);

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  std::uint32_t particle_count = 0;
  double density_value = -1.0;
  double temperature = -1.0;
  double vx = -1.0;
  double vy = -1.0;
  double vz = -1.0;
  double energy = -1.0;

  input.read(reinterpret_cast<char*>(&x), sizeof(x));
  input.read(reinterpret_cast<char*>(&y), sizeof(y));
  input.read(reinterpret_cast<char*>(&z), sizeof(z));
  input.read(reinterpret_cast<char*>(&particle_count), sizeof(particle_count));
  input.read(reinterpret_cast<char*>(&density_value), sizeof(density_value));
  input.read(reinterpret_cast<char*>(&temperature), sizeof(temperature));
  input.read(reinterpret_cast<char*>(&vx), sizeof(vx));
  input.read(reinterpret_cast<char*>(&vy), sizeof(vy));
  input.read(reinterpret_cast<char*>(&vz), sizeof(vz));
  input.read(reinterpret_cast<char*>(&energy), sizeof(energy));

  EXPECT_NEAR(x, 0.5, 1e-9);
  EXPECT_NEAR(y, 0.5, 1e-9);
  EXPECT_NEAR(z, 0.5, 1e-9);
  EXPECT_EQ(particle_count, 0U);
  EXPECT_NEAR(density_value, 0.0, 1e-12);
  EXPECT_NEAR(temperature, 0.0, 1e-12);
  EXPECT_NEAR(vx, 0.0, 1e-12);
  EXPECT_NEAR(vy, 0.0, 1e-12);
  EXPECT_NEAR(vz, 0.0, 1e-12);
  EXPECT_NEAR(energy, 0.0, 1e-12);

  char extra;
  EXPECT_FALSE(input.read(&extra, 1));
  EXPECT_TRUE(input.eof());
}

TEST(CellClusterIntegrationTest, RunsComputeAndInvokesBoundary) {
  const auto temp_dir =
      std::filesystem::temp_directory_path() / "mc3d_cluster_integration_test";
  ScopedWorkingDirectory guard(temp_dir);

  auto geometry = std::make_unique<mc3d::Geometry>();

  mc3d::CellCluster cluster;
  cluster.SetApex(mc3d::Point(0.0, 0.0, 0.0));
  cluster.SetSize(1.0, 1.0, 1.0);
  cluster.SetEndTime(0.05);
  cluster.SetSnapshotInterval(1.0);
  cluster.SetBinaryOutput(false);

  ASSERT_TRUE(cluster.Initialize(1, 1, 1, 10.0, kKn, kCu, std::move(geometry),
                                 kS, kAlpha, kTemperature));

  auto boundary = std::make_unique<mc3d::test::CountingBoundary>();
  auto* boundary_ptr = boundary.get();
  cluster.SetBoundaryCondition(std::move(boundary));

  cluster.Compute();

  EXPECT_NEAR(mc3d::CellClusterTestAccess::GetTime(cluster), 0.05, 1e-6);
  ASSERT_NE(boundary_ptr, nullptr);
  EXPECT_GT(boundary_ptr->call_count, 0);
  EXPECT_GT(boundary_ptr->last_dt, 0.0);
}
