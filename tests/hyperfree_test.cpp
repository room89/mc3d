#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>

#include "cell.h"
#include "giper_free_boundary.h"
#include "utils/utils.hpp"

TEST(HyperFreeSamplingTest, FluxWeightedNormalMomentsMatchMaxwellian) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  auto& random = *random_stream;
  constexpr double temperature = 1.4;
  constexpr int samples = 50000;
  const double sigma = std::sqrt(temperature);
  for (double drift : {-5., -1., 0., 1., 7.}) {
    SCOPED_TRACE(drift);
    random.Seed(12345);
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
      const double w = utils::IncomingNormalSpeed(drift, temperature, random);
      ASSERT_GT(w, 0);
      sum += w;
      squared += w * w;
    }
    EXPECT_NEAR(sum / samples, m2 / m1, 0.025 * sigma);
    EXPECT_NEAR(squared / samples, m3 / m1, 0.025 * (1 + m3 / m1));
  }
}

TEST(HyperFreeSamplingTest, FluxRetainsTheStrongOutflowTail) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  auto& random = *random_stream;
  // Integral of w*N(-10,1) over w>0, evaluated independently at high precision.
  EXPECT_NEAR(utils::IncomingFlux(-10, 1), 7.47456025458933e-25, 1e-35);
  EXPECT_NEAR(utils::IncomingFlux(10, 1), 10, 1e-12);
  EXPECT_THROW(utils::IncomingNormalSpeed(0, 0, random), std::invalid_argument);
}

TEST(HyperFreeSamplingTest, EveryFaceUsesIncomingVelocitiesAndUniformEntryTimes) {
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  auto& random = *random_stream;
  for (const mc3d::Point normal : {mc3d::Point(-1,0,0), {1,0,0}, {0,-1,0},
                                   {0,1,0}, {0,0,-1}, {0,0,1}}) {
    random.Seed(317);
    mc3d::Cell cell(random_stream);
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
  auto random_stream = std::make_shared<utils::RandomStream>(123);
  auto& random = *random_stream;
  random.Seed(19);
  auto body = std::make_unique<mc3d::Geometry>();
  body->CreateCube(0.02, 10, 0.1, 10);
  mc3d::Cell cell(random_stream);
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

TEST(HyperFreeSamplingTest, ReservoirDensityScalesFractionalFlux) {
  for(double ratio:{.37,1.,3.7}) for(double ux:{-1.,2.}) {
    std::deque<mc3d::Cell> cells;
    cells.emplace_back(std::make_shared<utils::RandomStream>(123));
    cells.front().SetApex({0,0,0});cells.front().SetSize(1,1,1);
    mc3d::HyperFreeBoundary boundary({0,0,0},{-1,0,0},20,0,1);
    boundary.SetReservoir(ratio,2,{ux,-3,1});boundary.BindCells(cells);
    size_t count=0;
    for(int i=0;i<1000;++i) count+=boundary.GenerateInflow(.001).transit.size();
    double expected=20*ratio*utils::IncomingFlux(ux,2);
    EXPECT_LE(double(count),expected+1e-10);
    EXPECT_GT(double(count),expected-1-1e-10);
    EXPECT_NEAR(boundary.TransportTimeLimit(),1/(3+2.),1e-14);
  }
}

TEST(HyperFreeSamplingTest, StochasticCountsHaveExpectedMeanAndAreReproducible) {
  auto sample = [](unsigned seed) {
    std::deque<mc3d::Cell> cells;
    cells.emplace_back(std::make_shared<utils::RandomStream>(seed));
    cells.front().SetApex({0,0,0}); cells.front().SetSize(1,1,1);
    mc3d::HyperFreeBoundary boundary({0,0,0},{-1,0,0},20,0,1);
    boundary.SetReservoir(1,1,{0,0,0}); boundary.BindCells(cells);
    boundary.SetStochasticRounding(true);
    std::vector<unsigned> counts;
    double total = 0, expected = 0, variance = 0;
    for (int i=0;i<20000;++i) {
      const double target = i%2 ? 2.7 : .2;
      const double dt = target/(20*utils::IncomingFlux(0,1));
      const auto n = boundary.GenerateInflow(dt).transit.size();
      EXPECT_TRUE(n==static_cast<unsigned>(target) || n==static_cast<unsigned>(target)+1);
      counts.push_back(static_cast<unsigned>(n));
      total += n; expected += target;
      const double fraction = target-std::floor(target);
      variance += fraction*(1-fraction);
    }
    EXPECT_NEAR(total, expected, 6*std::sqrt(variance));
    EXPECT_TRUE(boundary.GenerateInflow(0).transit.empty());
    return counts;
  };
  EXPECT_EQ(sample(917),sample(917));
  EXPECT_NE(sample(917),sample(918));
}

TEST(HyperFreeSamplingTest, ShockReservoirMomentumAndEnergyFluxes) {
  // Actual M=3.8 engineering-shock reservoirs; normalized mass and k_B are one.
  for (bool downstream : {false, true}) {
    SCOPED_TRACE(downstream);
    const double temperature = downstream ? 5.374515235457063 : 1.;
    const double density = downstream ? 3.3119266055045875 : 1.;
    const double ux = downstream ? 1.4812462622337137 : 4.905778905196061;
    const mc3d::Point normal(downstream ? 1 : -1, 0, 0);
    const double drift = downstream ? -ux : ux;
    const double sigma = std::sqrt(temperature);
    // Half-line Gaussian moments M_k = integral_0^inf w^k N(drift,T) dw.
    // Incoming velocity moments are M_(k+1)/M_1, not Gaussian volume moments.
    double m[6];
    m[0] = .5 * std::erfc(-drift / (std::sqrt(2.) * sigma));
    m[1] = drift * m[0] + sigma / std::sqrt(2 * std::numbers::pi) *
        std::exp(-drift * drift / (2 * temperature));
    for (int k = 2; k <= 5; ++k)
      m[k] = drift * m[k-1] + (k-1) * temperature * m[k-2];
    const double expected_w = m[2]/m[1];
    const double variance_w = m[3]/m[1] - expected_w*expected_w;
    const double expected_energy = .5 * (m[3]/m[1] + 2*temperature);
    const double variance_energy = .25 * (m[5]/m[1] -
        std::pow(m[3]/m[1],2) + 4*temperature*temperature);
    auto stream = std::make_shared<utils::RandomStream>(20261116);
    std::deque<mc3d::Cell> cells;
    cells.emplace_back(stream);
    cells.front().SetApex({-.5,-.5,-.5});
    cells.front().SetSize(1,1,1);
    mc3d::HyperFreeBoundary boundary(normal*.5, normal, 160, 0, temperature);
    boundary.SetReservoir(density,temperature,{ux,0,0});
    boundary.BindCells(cells);
    double expected_count = 0, sum_w = 0, sum_energy = 0;
    size_t count = 0;
    for (int step = 0; step < 1000; ++step) {
      // Unequal steps exercise the fractional carry, not only the sampler.
      const double target = step % 2 ? 173.375 : 226.125;
      const double dt = target / (160*density*m[1]);
      auto arrivals = boundary.GenerateInflow(dt).transit;
      expected_count += target;
      count += arrivals.size();
      for (const auto& p : arrivals) {
        const double w = -(p.velocity*normal);
        ASSERT_GT(w,0);
        sum_w += w;
        sum_energy += .5*(p.velocity*p.velocity);
      }
    }
    EXPECT_LE(std::abs(double(count)-expected_count),1.00000001);
    // Six analytic standard errors; fixed seed, no tuned tolerance or resampling.
    EXPECT_NEAR(sum_w/count,expected_w,6*std::sqrt(variance_w/count));
    EXPECT_NEAR(sum_energy/count,expected_energy,6*std::sqrt(variance_energy/count));
  }
}
