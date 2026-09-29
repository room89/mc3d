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
    boundary->AddCell(cells);
    return boundary;
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
    boundary->BoundaryCondition(outgoing, DurationForFlux(0.2025));
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
      boundary->BoundaryCondition(outgoing, DurationForFlux(flux));
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
    boundary->BoundaryCondition(outgoing, DurationForFlux(flux));
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
    boundary->BoundaryCondition(outgoing, DurationForFlux(0.1717));
    a->BoundaryCondition(arrivals_a, DurationForFlux(0.1717));
    b->BoundaryCondition(arrivals_b, DurationForFlux(0.1717));
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
  coarse->BoundaryCondition(coarse_arrivals, DurationForFlux(30.5));
  for (int step = 0; step < 100; ++step) {
    fine->BoundaryCondition(fine_arrivals, DurationForFlux(0.305));
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
