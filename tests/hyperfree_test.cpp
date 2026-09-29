#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>

#include "cell.h"
#include "giper_free_boundary.h"
#include "utils/utils.hpp"

TEST(HyperFreeSamplingTest, FluxWeightedNormalMomentsMatchMaxwellian) {
  constexpr double temperature = 1.4;
  constexpr int samples = 50000;
  const double sigma = std::sqrt(temperature);
  for (double drift : {-5., -1., 0., 1., 7.}) {
    SCOPED_TRACE(drift);
    utils::RandomEngine().seed(12345);
    const double phi = std::exp(-drift * drift / (2 * temperature)) /
                       std::sqrt(2 * std::numbers::pi);
    const double cdf = 0.5 * std::erfc(-drift / (std::sqrt(2.) * sigma));
    const double m1 = drift * cdf + sigma * phi;
    const double m2 = (drift * drift + temperature) * cdf + drift * sigma * phi;
    const double m3 = (drift * drift * drift + 3 * drift * temperature) * cdf +
                     (drift * drift + 2 * temperature) * sigma * phi;
    EXPECT_NEAR(utils::IncomingFlux(drift, temperature), m1, 1e-14);
    double sum = 0, squared = 0;
    for (int i = 0; i < samples; ++i) {
      const double w = utils::IncomingNormalSpeed(drift, temperature);
      ASSERT_GT(w, 0);
      sum += w;
      squared += w * w;
    }
    EXPECT_NEAR(sum / samples, m2 / m1, 0.025 * sigma);
    EXPECT_NEAR(squared / samples, m3 / m1, 0.025 * (1 + m3 / m1));
  }
}

TEST(HyperFreeSamplingTest, FluxRetainsTheStrongOutflowTail) {
  // Integral of w*N(-10,1) over w>0, evaluated independently at high precision.
  EXPECT_NEAR(utils::IncomingFlux(-10, 1), 7.47456025458933e-25, 1e-35);
  EXPECT_NEAR(utils::IncomingFlux(10, 1), 10, 1e-12);
  EXPECT_THROW(utils::IncomingNormalSpeed(0, 0), std::invalid_argument);
}

TEST(HyperFreeSamplingTest, EveryFaceUsesIncomingVelocitiesAndUniformEntryTimes) {
  for (const mc3d::Point normal : {mc3d::Point(-1,0,0), {1,0,0}, {0,-1,0},
                                   {0,1,0}, {0,0,-1}, {0,0,1}}) {
    utils::RandomEngine().seed(317);
    mc3d::Cell cell;
    cell.SetApex({-0.5,-0.5,-0.5});
    cell.SetSize(1,1,1);
    std::vector<mc3d::Particle> arrivals;
    constexpr double dt = 0.2;
    const mc3d::Point drift(3,-1,0.5);
    cell.GenerateHyperFreeRandom(10000, drift, 1, normal, dt, arrivals);
    ASSERT_EQ(arrivals.size(), 10000U);
    EXPECT_EQ(cell.GetParticleCount(), 0U);
    double mean_age = 0;
    mc3d::Point mean_velocity(0,0,0);
    for (const auto& p : arrivals) {
      ASSERT_LT(p.velocity * normal, 0);
      const double age = ((p.position - normal * 0.5) * normal) / (p.velocity * normal);
      EXPECT_GE(age, -1e-14);
      EXPECT_LE(age, dt + 1e-14);
      const auto source = p.position - p.velocity * age;
      EXPECT_LE(std::abs(source.x), 0.5 + 1e-12);
      EXPECT_LE(std::abs(source.y), 0.5 + 1e-12);
      EXPECT_LE(std::abs(source.z), 0.5 + 1e-12);
      mean_age += age;
      mean_velocity += p.velocity;
    }
    EXPECT_NEAR(mean_age / arrivals.size(), dt / 2, 0.003);
    const auto tangent_error = mean_velocity / double(arrivals.size()) - drift;
    const auto projected = tangent_error - normal * (tangent_error * normal);
    EXPECT_LT(projected.Mod(), 0.05);
    // In particular Y/Z planes through the origin must not be mistaken for X.
    mc3d::HyperFreeBoundary boundary({0,0,0}, normal, 10, 2, 1);
    EXPECT_DOUBLE_EQ(boundary.GetNormal() * normal, 1);
  }
}

TEST(HyperFreeSamplingTest, NewArrivalsCannotTunnelThroughNearbyBody) {
  utils::RandomEngine().seed(19);
  auto body = std::make_unique<mc3d::Geometry>();
  body->CreateCube(0.02, 10, 0.1, 10);
  mc3d::Cell cell;
  cell.SetApex({0,-0.5,-0.5});
  cell.SetSize(1,1,1);
  cell.SetParameters(0,0,1);
  ASSERT_TRUE(cell.Initialize(0, body));
  std::vector<mc3d::Particle> arrivals;
  cell.GenerateHyperFreeRandom(1000, {7,0,0}, 0.1, {-1,0,0}, 0.1, arrivals);
  int reflected = 0;
  for (const auto& particle : arrivals) {
    EXPECT_FALSE(body->IsInnerPoint(particle.position));
    EXPECT_LT(particle.position.x, 0.02);
    reflected += particle.velocity.x < 0;
  }
  EXPECT_GT(reflected, 900);
}
