#include <gtest/gtest.h>

#include <cmath>
#include <deque>
#include <memory>
#include <numbers>
#include <vector>

#include "free_boundary.h"
#include "giper_free_boundary.h"

namespace {

class OpenBoundaryInflowTest : public testing::TestWithParam<bool> {
 protected:
  static constexpr unsigned kParticlesPerCell = 60;
  static constexpr double kTemperature = 1.2;

  std::unique_ptr<mc3d::FreeBoundary> MakeBoundary(std::deque<mc3d::Cell>& cells) {
    std::unique_ptr<mc3d::FreeBoundary> boundary;
    if (GetParam()) {
      boundary = std::make_unique<mc3d::HyperFreeBoundary>(
          mc3d::Point(0, 0, -0.5), mc3d::Point(0, 0, -1),
          kParticlesPerCell, 10, kTemperature);
    } else {
      boundary = std::make_unique<mc3d::FreeBoundary>(
          mc3d::Point(0, 0, -0.5), mc3d::Point(0, 0, -1),
          kParticlesPerCell, 10, kTemperature);
    }
    boundary->BindCells(cells);
    return boundary;
  }

  static void ApplyInflow(mc3d::Boundary& boundary, std::vector<mc3d::Particle>& arrivals, double dt) {
    auto result = boundary.GenerateInflow(dt);
    for (auto& batch : result.local)
      batch.cell.get().AppendGeneratedParticles(std::move(batch.particles));
    arrivals.insert(arrivals.end(), result.transit.begin(), result.transit.end());
  }

  static void AddCell(std::deque<mc3d::Cell>& cells, double width = 1) {
    cells.emplace_back();
    cells.back().SetApex({double(cells.size() - 1), 0, -0.5});
    cells.back().SetSize(1, 1, width);
  }

  // Zero normal drift: the Maxwellian incoming flux is n * sqrt(T / 2pi).
  static double DurationForFlux(double count) {
    return count /
           (kParticlesPerCell * std::sqrt(kTemperature / (2 * std::numbers::pi)));
  }

  void ExpectBudget(const mc3d::Cell& cell, double expected,
                    const std::vector<mc3d::Particle>& arrivals) {
    const double actual = cell.GetParticleCount() + arrivals.size();
    EXPECT_LE(actual, expected + 1e-10);
    EXPECT_LT(expected - actual, (GetParam() ? 1 : 3) + 1e-10);
  }
};

TEST_P(OpenBoundaryInflowTest, SubParticleFluxIsNotLostAtSmallTimeSteps) {
  std::deque<mc3d::Cell> cells;
  AddCell(cells);
  auto boundary = MakeBoundary(cells);
  std::vector<mc3d::Particle> outgoing;
  for (int step = 1; step <= 100; ++step) {
    ApplyInflow(*boundary, outgoing, DurationForFlux(0.2025));
    ExpectBudget(cells.front(), step * 0.2025, outgoing);
  }
  EXPECT_GT(cells.front().GetParticleCount() + outgoing.size(), 17U);
}

TEST_P(OpenBoundaryInflowTest, OneAndTwoParticleInflowsAreNotDiscarded) {
  for (const double flux : {1.25, 2.25}) {
    SCOPED_TRACE(flux);
    std::deque<mc3d::Cell> cells;
    AddCell(cells);
    auto boundary = MakeBoundary(cells);
    std::vector<mc3d::Particle> outgoing;
    for (int step = 1; step <= 10; ++step) {
      ApplyInflow(*boundary, outgoing, DurationForFlux(flux));
      ExpectBudget(cells.front(), step * flux, outgoing);
    }
  }
}

TEST_P(OpenBoundaryInflowTest, FractionalFluxSurvivesVariableTimeSteps) {
  std::deque<mc3d::Cell> cells;
  AddCell(cells);
  auto boundary = MakeBoundary(cells);
  std::vector<mc3d::Particle> outgoing;
  double expected = 0;
  for (const double flux : {0.25, 3.2, 1.1, 0.0, 2.3, 4.7, 0.1}) {
    expected += flux;
    ApplyInflow(*boundary, outgoing, DurationForFlux(flux));
    ExpectBudget(cells.front(), expected, outgoing);
  }
}

TEST_P(OpenBoundaryInflowTest, EachCellAccumulatesItsOwnFlux) {
  std::deque<mc3d::Cell> cells;
  AddCell(cells, 1);
  AddCell(cells, 0.5);
  auto boundary = MakeBoundary(cells);
  std::deque<mc3d::Cell> reference_a, reference_b;
  AddCell(reference_a, 1);
  AddCell(reference_b, 0.5);
  auto a = MakeBoundary(reference_a), b = MakeBoundary(reference_b);
  std::vector<mc3d::Particle> outgoing;
  std::vector<mc3d::Particle> arrivals_a, arrivals_b;
  for (int step = 1; step <= 100; ++step) {
    ApplyInflow(*boundary, outgoing, DurationForFlux(0.1717));
    ApplyInflow(*a, arrivals_a, DurationForFlux(0.1717));
    ApplyInflow(*b, arrivals_b, DurationForFlux(0.1717));
    ExpectBudget(reference_a[0], step * 0.1717, arrivals_a);
    ExpectBudget(reference_b[0], step * 0.3434, arrivals_b);
    EXPECT_EQ(cells[0].GetParticleCount() + cells[1].GetParticleCount() + outgoing.size(),
              reference_a[0].GetParticleCount() + reference_b[0].GetParticleCount() +
                  arrivals_a.size() + arrivals_b.size());
  }
}

TEST_P(OpenBoundaryInflowTest, SubdividingTimeKeepsIntegratedInflow) {
  std::deque<mc3d::Cell> coarse_cells, fine_cells;
  AddCell(coarse_cells);
  AddCell(fine_cells);
  auto coarse = MakeBoundary(coarse_cells);
  auto fine = MakeBoundary(fine_cells);
  std::vector<mc3d::Particle> coarse_arrivals, fine_arrivals;
  ApplyInflow(*coarse, coarse_arrivals, DurationForFlux(30.5));
  for (int step = 0; step < 100; ++step) {
    ApplyInflow(*fine, fine_arrivals, DurationForFlux(0.305));
  }
  EXPECT_EQ(coarse_cells.front().GetParticleCount() + coarse_arrivals.size(),
            fine_cells.front().GetParticleCount() + fine_arrivals.size());
  EXPECT_EQ(fine_cells.front().GetParticleCount() + fine_arrivals.size(), 30U);
}

INSTANTIATE_TEST_SUITE_P(
    FreeAndHyperFree, OpenBoundaryInflowTest, testing::Bool(),
    [](const testing::TestParamInfo<bool>& info) {
      return info.param ? "HyperFree" : "Free";
    });

}  // namespace

TEST(BoundaryPhasesTest, InflowReturnsBatchesWithoutChangingCells) {
  std::deque<mc3d::Cell> cells;
  cells.emplace_back(std::make_shared<utils::RandomStream>(7));
  cells.front().SetApex({0,0,0});
  cells.front().SetSize(1,1,1);
  mc3d::FreeBoundary free({0,0,0}, {-1,0,0}, 60, 0, 1);
  mc3d::HyperFreeBoundary hyper({0,0,0}, {-1,0,0}, 60, 0, 1);
  free.BindCells(cells);
  hyper.BindCells(cells);
  const double dt = 3.25 / (60 * utils::IncomingFlux(0, 1));
  const auto local = free.GenerateInflow(dt);
  ASSERT_EQ(local.local.size(), 1u);
  EXPECT_EQ(&local.local.front().cell.get(), &cells.front());
  EXPECT_EQ(local.local.front().particles.size(), 3u);
  EXPECT_TRUE(local.transit.empty());
  const auto transit = hyper.GenerateInflow(dt);
  EXPECT_TRUE(transit.local.empty());
  EXPECT_EQ(transit.transit.size(), 3u);
  EXPECT_EQ(cells.front().GetParticleCount(), 0u);
  std::vector<mc3d::Particle> outgoing{{{2,0,0},{1,0,0}}};
  free.ProcessOutgoing(outgoing, dt);
  EXPECT_EQ(outgoing.size(), 1u);
  EXPECT_DOUBLE_EQ(outgoing.front().position.x, 2);
}
