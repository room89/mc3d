#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <thread>
#include <vector>

#include "geometry.h"
#include "inner_boundary.h"
#include "simulation_runner.h"

namespace {
class ThermalBody : public mc3d::Geometry {
 public:
  mc3d::Polygon& Face(std::size_t i) { return *poligons.at(i); }
};

double Energy(const std::vector<mc3d::Particle>& particles) {
  double energy = 0;
  for (const auto& p : particles) energy += 0.5 * (p.velocity * p.velocity);
  return energy;
}

TEST(SurfaceHeatingTest, ReflectionsConserveGasPlusSurfaceEnergyAndMomentum) {
  for (double initial_temperature : {1., 100.}) {
    ThermalBody body;
    body.CreateCube(-0.5, 1, 1, 1);
    body.ConfigureSurfaceHeating(initial_temperature, 100, 0.001);
    mc3d::InnerBoundary boundary;
    boundary.SetGeometry(body);
    // Local polygon temperature overrides the fixed-wall fallback.
    boundary.SetWallTemperature(1000);
    utils::RandomStream random(123);
    std::vector<mc3d::Particle> particles(2000, {{-0.6,0.1,0.2},{10,0,0}});
    const double incoming = Energy(particles);
    boundary.BoundaryCondition(particles, 0.011, random);
    const double gas_change = 0.001 * (Energy(particles) - incoming);
    mc3d::Point gas_momentum_change(0,0,0);
    for (const auto& p : particles)
      gas_momentum_change += p.velocity - mc3d::Point(10,0,0);
    const auto impulse = body.TakeStepImpulse();
    EXPECT_NEAR((impulse + gas_momentum_change).Mod(), 0, 1e-8);
    EXPECT_DOUBLE_EQ(body.TakeStepImpulse().Mod(), 0);
    body.UpdateSurfaceTemperatures();
    double surface_change = 0, recorded_heat = 0;
    int changed = 0;
    for (std::size_t i = 0; i < body.PolygonCount(); ++i) {
      const auto& p = body.GetPolygon(i);
      surface_change += p.heat_capacity * (p.wall_temperature - initial_temperature);
      recorded_heat += p.absorbed_heat;
      if (p.wall_temperature != initial_temperature) ++changed;
    }
    EXPECT_EQ(changed, 1);
    EXPECT_NEAR(surface_change + gas_change, 0, 1e-9);
    EXPECT_NEAR(recorded_heat, surface_change, 1e-9);
    if (initial_temperature == 1) EXPECT_GT(surface_change, 0);
    else EXPECT_LT(surface_change, 0);
    body.UpdateSurfaceTemperatures();
    double heat_after_empty_step = 0;
    for (std::size_t i = 0; i < body.PolygonCount(); ++i)
      heat_after_empty_step += body.GetPolygon(i).absorbed_heat;
    EXPECT_DOUBLE_EQ(heat_after_empty_step, recorded_heat);
  }
}

TEST(SurfaceHeatingTest, CapacityFollowsAreaAndParticleWeightNormalizesHeating) {
  for (int copies : {1, 10}) {
    ThermalBody body;
    body.CreateCube(-0.5, 1, 2, 3);
    body.ConfigureSurfaceHeating(1, 44, 0.1 / copies);
    double capacity = 0;
    for (std::size_t i = 0; i < body.PolygonCount(); ++i) {
      const auto& p = body.GetPolygon(i);
      const double area = 0.5 * (p.p2-p.p1).Cross(p.p3-p.p1).Mod();
      EXPECT_NEAR(p.heat_capacity / area, 2, 1e-12);
      capacity += p.heat_capacity;
    }
    EXPECT_NEAR(capacity, 44, 1e-12);
    for (int n = 0; n < copies; ++n)
      body.AccumulateForce(body.Face(0), {0,0,0}, 20);
    body.UpdateSurfaceTemperatures();
    EXPECT_NEAR(body.Face(0).absorbed_heat, 2, 1e-12);
    EXPECT_NEAR(body.Face(0).wall_temperature, 1+2/body.Face(0).heat_capacity, 1e-12);
  }
}

TEST(SurfaceHeatingTest, RejectsInvalidSetupAndDoesNotClampExhaustedEnergy) {
  ThermalBody body;
  EXPECT_THROW(body.ConfigureSurfaceHeating(1,1,1), std::invalid_argument);
  body.CreateCube();
  for (double invalid : {0., -1., std::numeric_limits<double>::infinity(),
                         std::numeric_limits<double>::quiet_NaN()}) {
    EXPECT_THROW(body.ConfigureSurfaceHeating(1,invalid,1), std::invalid_argument);
    EXPECT_THROW(body.ConfigureSurfaceHeating(invalid,1,1), std::invalid_argument);
    EXPECT_THROW(body.ConfigureSurfaceHeating(1,1,invalid), std::invalid_argument);
  }
  body.ConfigureSurfaceHeating(1,12,1);
  body.AccumulateForce(body.Face(0), {0,0,0}, 1);
  body.AccumulateForce(body.Face(1), {0,0,0}, -100);
  EXPECT_THROW(body.UpdateSurfaceTemperatures(), std::runtime_error);
  EXPECT_DOUBLE_EQ(body.Face(0).wall_temperature, 1);
  EXPECT_DOUBLE_EQ(body.Face(1).wall_temperature, 1);
  EXPECT_DOUBLE_EQ(body.Face(1).pending_heat, -100);
}

TEST(SurfaceHeatingTest, ConcurrentImpactsPreserveHeatAndImpulse) {
  ThermalBody body;
  body.CreateCube();
  body.ConfigureSurfaceHeating(1,12,0.01);
  std::vector<std::thread> threads;
  for (int i = 0; i < 4; ++i)
    threads.emplace_back([&] {
      for (int n = 0; n < 1000; ++n)
        body.AccumulateForce(body.Face(0), {2,-1,3}, 4);
    });
  for (auto& thread : threads) thread.join();
  const auto impulse = body.TakeStepImpulse();
  EXPECT_DOUBLE_EQ(impulse.x, 8000);
  EXPECT_DOUBLE_EQ(impulse.y, -4000);
  EXPECT_DOUBLE_EQ(impulse.z, 12000);
  body.UpdateSurfaceTemperatures();
  EXPECT_NEAR(body.Face(0).absorbed_heat, 160, 1e-9);
}

TEST(SurfaceHeatingTest, FixedWallStillTracksMomentumWithoutHeating) {
  ThermalBody body;
  body.CreateCube();
  body.AccumulateForce(body.Face(0), {2,3,4}, 100);
  body.UpdateSurfaceTemperatures();
  EXPECT_FALSE(body.SurfaceHeatingEnabled());
  EXPECT_DOUBLE_EQ(body.Face(0).pending_heat, 0);
  EXPECT_DOUBLE_EQ(body.Face(0).wall_temperature, 1);
  EXPECT_DOUBLE_EQ(body.TakeStepImpulse().x, 2);
}

TEST(BodyLoadsTest, RejectsUndefinedCoefficientNormalization) {
  mc3d::CellCluster cluster(1, 1);
  mc3d::SimulationRunner runner(cluster);
  EXPECT_THROW(runner.SetForceReference(0,1,0), std::invalid_argument);
  EXPECT_THROW(runner.SetForceReference(1,0,0), std::invalid_argument);
  EXPECT_THROW(runner.SetForceReference(1,1,std::numeric_limits<double>::infinity()), std::invalid_argument);
  EXPECT_NO_THROW(runner.SetForceReference(0.2,7,0.1));
}
}  // namespace
