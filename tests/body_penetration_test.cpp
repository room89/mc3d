#include <gtest/gtest.h>

#include <array>
#include <memory>
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
