#include <gtest/gtest.h>

#include <array>
#include <memory>
#include <thread>
#include <vector>

#include "cell.h"
#include "geometry.h"
#include "inner_boundary.h"
#include "particle.h"
#include "poligon.h"

namespace {

class MeshGeometry : public mc3d::Geometry {
 public:
  void Triangle(mc3d::Point a, mc3d::Point b, mc3d::Point c,
                mc3d::Point normal) {
    auto polygon = std::make_unique<mc3d::Polygon>(a, b, c, normal);
    polygon->force = {0, 0, 0};
    poligons.push_back(std::move(polygon));
  }
};

std::unique_ptr<MeshGeometry> MakeConcavePrism() {
  auto geometry = std::make_unique<MeshGeometry>();
  // Extrude the L-shaped polygon (0,0)-(2,0)-(2,1)-(1,1)-(1,2)-(0,2).
  const std::array<mc3d::Point, 6> outline{
      {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}}};
  for (std::size_t i = 0; i < outline.size(); ++i) {
    auto a = outline[i];
    auto b = outline[(i + 1) % outline.size()];
    auto top_a = a + mc3d::Point(0, 0, 1);
    auto top_b = b + mc3d::Point(0, 0, 1);
    mc3d::Point normal(b.y - a.y, a.x - b.x, 0);
    geometry->Triangle(a, b, top_b, normal);
    geometry->Triangle(a, top_b, top_a, normal);
  }
  const std::array<std::array<mc3d::Point, 4>, 3> squares{
      {{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}},
       {{{1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}}},
       {{{0, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}}}}};
  for (const auto& square : squares) {
    geometry->Triangle(square[0], square[2], square[1], {0, 0, -1});
    geometry->Triangle(square[0], square[3], square[2], {0, 0, -1});
    const auto rise = mc3d::Point(0, 0, 1);
    geometry->Triangle(square[0] + rise, square[1] + rise, square[2] + rise,
                       {0, 0, 1});
    geometry->Triangle(square[0] + rise, square[2] + rise, square[3] + rise,
                       {0, 0, 1});
  }
  return geometry;
}

}  // namespace

TEST(BodyGeometryTest, ClassifiesConcaveClosedMesh) {
  auto body = MakeConcavePrism();
  EXPECT_TRUE(body->IsInnerPoint({0.5, 0.5, 0.5}));
  EXPECT_TRUE(body->IsInnerPoint({1.5, 0.5, 0.5}));
  EXPECT_TRUE(body->IsInnerPoint({0.5, 1.5, 0.5}));
  EXPECT_FALSE(body->IsInnerPoint({1.5, 1.5, 0.5}));
  EXPECT_FALSE(body->IsInnerPoint({-0.5, 0.5, 0.5}));
  EXPECT_FALSE(body->IsInnerPoint({1, 1.5, 0.5}));
}

TEST(BodyGeometryTest, RecoversStartingParticleFromDeepInside) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  ASSERT_TRUE(body.IsInnerPoint({0, 0, 0}));
  const auto exterior = body.ExteriorPoint({0, 0, 0});
  ASSERT_TRUE(exterior.has_value());
  EXPECT_FALSE(body.IsInnerPoint(*exterior));
}

TEST(BodyCollisionTest, HighSpeedCrossingCannotTunnelThroughCube) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles{{{-1, 0, 0}, {4, 0, 0}},
                                        {{-1, 0.49, 0.49}, {4, 0, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.5), 0);
  for (const auto& particle : particles) {
    EXPECT_FALSE(body.IsInnerPoint(particle.position));
    EXPECT_LT(particle.position.x, -0.5);
  }
}

TEST(BodyCollisionTest, ReflectsAtConcaveNotchWithoutPenetration) {
  auto body = MakeConcavePrism();
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(*body);
  std::vector<mc3d::Particle> particles{{{1.5, 1.5, 0.5}, {-3, 0, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.3), 0);
  EXPECT_FALSE(body->IsInnerPoint(particles.front().position));
  EXPECT_GT(particles.front().velocity.x, 0);
}

TEST(BodyCollisionTest, AdjacentCellDetectsSurfaceOutsideLocalPolygonList) {
  auto body = std::make_unique<mc3d::Geometry>();
  body->CreateCube(-0.5, 1, 1, 1);
  mc3d::Cell cell;
  cell.SetApex({-1.5, -0.5, -0.5});
  cell.SetSize(1, 1, 1);
  cell.SetParameters(0, 0, 1);
  ASSERT_TRUE(cell.Initialize(0, body));
  cell.SetDt(0.5);
  mc3d::Particle particle({-1, 0, 0}, {2, 0, 0});
  ASSERT_TRUE(cell.TryAcceptParticle(particle));
  cell.Calculate();
  EXPECT_EQ(cell.GetParticleCount(), 1U);
  EXPECT_EQ(cell.CountInnerParticles(*body), 0U);
  EXPECT_LT(cell.GetParticleMassCenter().x, -0.5);
}

TEST(BodyCollisionTest, RepeatedCrossingsKeepEveryParticleOutside) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles;
  for (int i = 0; i < 60; ++i) {
    const double offset = 0.4 * ((i % 9) - 4) / 4.0;
    const double speed = 1.0 + (i % 7);
    switch (i % 6) {
      case 0:
        particles.emplace_back(mc3d::Point(-1, offset, 0),
                               mc3d::Point(speed, 0, 0));
        break;
      case 1:
        particles.emplace_back(mc3d::Point(1, offset, 0),
                               mc3d::Point(-speed, 0, 0));
        break;
      case 2:
        particles.emplace_back(mc3d::Point(offset, -1, 0),
                               mc3d::Point(0, speed, 0));
        break;
      case 3:
        particles.emplace_back(mc3d::Point(offset, 1, 0),
                               mc3d::Point(0, -speed, 0));
        break;
      case 4:
        particles.emplace_back(mc3d::Point(offset, 0, -1),
                               mc3d::Point(0, 0, speed));
        break;
      default:
        particles.emplace_back(mc3d::Point(offset, 0, 1),
                               mc3d::Point(0, 0, -speed));
        break;
    }
  }
  for (int step = 0; step < 20; ++step) {
    SCOPED_TRACE(step);
    ASSERT_EQ(boundary.BoundaryCondition(particles, 0.1), 0);
    ASSERT_EQ(particles.size(), 60U);
    for (std::size_t i = 0; i < particles.size(); ++i) {
      SCOPED_TRACE(i);
      EXPECT_FALSE(body.IsInnerPoint(particles[i].position));
    }
  }
}

TEST(BodyCollisionTest, ReversedMeshNormalsStillReflectIncomingParticle) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  body.ReverseNormals();
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles{{{-1, 0, 0}, {2, 0, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.3), 0);
  EXPECT_FALSE(body.IsInnerPoint(particles.front().position));
  EXPECT_LT(particles.front().velocity.x, 0);
}

TEST(BodyCollisionTest, SurfaceStartReflectsInsteadOfBeingRecoveredInside) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles{{{-0.5, 0, 0}, {1, 0, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.1), 0);
  EXPECT_FALSE(body.IsInnerPoint(particles.front().position));
  EXPECT_LT(particles.front().velocity.x, 0);
}

TEST(BodyGeometryTest, SegmentChoosesFirstFaceAndRejectsOutgoingStart) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  const auto hit = body.FirstIntersection({-1, 0, 0}, {2, 0, 0});
  ASSERT_TRUE(hit.has_value());
  EXPECT_NEAR(hit->fraction, 0.25, 1e-12);
  EXPECT_NEAR(hit->polygon->GetNormal().x, -1, 1e-12);
  EXPECT_FALSE(body.FirstIntersection({-0.5, 0, 0}, {-1, 0, 0}));
}

TEST(BodyGeometryTest, BoundsIncludeNegativeCoordinatesAndTranslatedBody) {
  mc3d::Geometry body;
  body.CreateCube(-3, 1, 1, 1);
  body.Move({-1, -2, -3});
  const auto [lower, upper] = body.Bounds();
  EXPECT_DOUBLE_EQ(lower.x, -4);
  EXPECT_DOUBLE_EQ(upper.x, -3);
  EXPECT_DOUBLE_EQ(lower.y, -2.5);
  EXPECT_DOUBLE_EQ(upper.y, -1.5);
  EXPECT_DOUBLE_EQ(lower.z, -3.5);
  EXPECT_DOUBLE_EQ(upper.z, -2.5);
}

TEST(BodyCollisionTest,
     StartingInsideIsRecoveredByBoundaryWithoutParticleLoss) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles{{{0, 0, 0}, {0, 0, 0}}};
  ASSERT_TRUE(body.IsInnerPoint(particles.front().position));
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.1), 0);
  ASSERT_EQ(particles.size(), 1U);
  EXPECT_FALSE(body.IsInnerPoint(particles.front().position));
  EXPECT_DOUBLE_EQ(particles.front().velocity.Mod(), 0);
}

TEST(BodyCollisionTest, DistantTrajectoryKeepsVelocityAndCreatesNoForce) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  std::vector<mc3d::Particle> particles{{{2, 2, 2}, {1, -0.25, 0.5}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.2), 0);
  EXPECT_NEAR(particles.front().position.x, 2.2, 1e-12);
  EXPECT_NEAR(particles.front().position.y, 1.95, 1e-12);
  EXPECT_NEAR(particles.front().position.z, 2.1, 1e-12);
  EXPECT_DOUBLE_EQ(particles.front().velocity.x, 1);
  EXPECT_DOUBLE_EQ(particles.front().velocity.y, -0.25);
  EXPECT_DOUBLE_EQ(particles.front().velocity.z, 0.5);
  for (std::size_t i = 0; i < body.PolygonCount(); ++i) {
    EXPECT_DOUBLE_EQ(body.GetPolygon(i).force.Mod(), 0);
  }
}

TEST(BodyCollisionTest, ConcurrentHitsAccumulateTheFullWallImpulse) {
  mc3d::Geometry body;
  body.CreateCube(-0.5, 1, 1, 1);
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(body);
  constexpr int kThreads = 8;
  constexpr int kParticlesPerThread = 100;
  std::array<mc3d::Point, kThreads> expected;
  std::array<bool, kThreads> stayed_outside{};
  std::vector<std::thread> workers;
  for (int thread_index = 0; thread_index < kThreads; ++thread_index) {
    workers.emplace_back([&, thread_index] {
      std::vector<mc3d::Particle> particles;
      for (int j = 0; j < kParticlesPerThread; ++j) {
        particles.emplace_back(mc3d::Point(-1, 0.1, 0.1), mc3d::Point(2, 0, 0));
      }
      boundary.BoundaryCondition(particles, 0.3);
      mc3d::Point impulse(0, 0, 0);
      bool outside = particles.size() == kParticlesPerThread;
      for (const auto& particle : particles) {
        impulse += mc3d::Point(2, 0, 0) - particle.velocity;
        outside &= !body.IsInnerPoint(particle.position);
      }
      expected[thread_index] = impulse;
      stayed_outside[thread_index] = outside;
    });
  }
  for (auto& worker : workers) worker.join();
  mc3d::Point total_expected(0, 0, 0);
  mc3d::Point total_force(0, 0, 0);
  for (int i = 0; i < kThreads; ++i) {
    EXPECT_TRUE(stayed_outside[i]);
    total_expected += expected[i];
  }
  for (std::size_t i = 0; i < body.PolygonCount(); ++i) {
    total_force += body.GetPolygon(i).force;
  }
  EXPECT_NEAR(total_force.x, total_expected.x, 1e-8);
  EXPECT_NEAR(total_force.y, total_expected.y, 1e-8);
  EXPECT_NEAR(total_force.z, total_expected.z, 1e-8);
}

TEST(BodyCollisionTest, RecoversFromConcaveInteriorTowardNotch) {
  auto body = MakeConcavePrism();
  const mc3d::Point start(0.9, 1.5, 0.5);
  ASSERT_TRUE(body->IsInnerPoint(start));
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(*body);
  std::vector<mc3d::Particle> particles{{start, {0, 0, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 0.2), 0);
  ASSERT_EQ(particles.size(), 1U);
  EXPECT_FALSE(body->IsInnerPoint(particles.front().position));
  EXPECT_GT(particles.front().position.x, 1);
}

TEST(BodyCollisionTest, ConcaveCornerHitDoesNotLeaveParticleInside) {
  auto body = MakeConcavePrism();
  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(*body);
  std::vector<mc3d::Particle> particles{{{1.5, 1.5, 0.5}, {-1, -1, 0}}};
  ASSERT_EQ(boundary.BoundaryCondition(particles, 1), 0);
  ASSERT_EQ(particles.size(), 1U);
  EXPECT_FALSE(body->IsInnerPoint(particles.front().position));
}
