#include <gtest/gtest.h>
#include <inner_boundary.h>
#include <particle.h>
#include <poligon.h>

#include <utils/utils.hpp>
#include <vector>

namespace mc3d::test {

class GeometryAccessor : public Geometry {
 public:
  using Geometry::Geometry;

  void AddPolygon(const Polygon& polygon) {
    poligons.push_back(std::make_unique<Polygon>(polygon));
  }

  Polygon& PolygonAt(size_t index) { return *poligons.at(index); }

  std::vector<Polygon*> RawPolygons() {
    std::vector<Polygon*> result;
    result.reserve(poligons.size());
    for (auto& polygon : poligons) {
      result.push_back(polygon.get());
    }
    return result;
  }
};

namespace {

Point MakePoint(double x, double y, double z) { return Point(x, y, z); }

constexpr double kTolerance = 1e-12;

}  // namespace

}  // namespace mc3d::test

namespace {

using mc3d::Particle;
using mc3d::Point;
using mc3d::test::kTolerance;

Point Difference(const Point& a, const Point& b) {
  return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}

}  // namespace

TEST(InnerBoundaryCollisionTest, ParticleCollidesWithPolygon) {
  mc3d::test::GeometryAccessor geometry;
  mc3d::Polygon polygon(mc3d::test::MakePoint(0.0, -0.5, -0.5),
                        mc3d::test::MakePoint(0.0, 0.5, -0.5),
                        mc3d::test::MakePoint(0.0, -0.5, 0.5),
                        mc3d::test::MakePoint(-1.0, 0.0, 0.0));
  polygon.force = mc3d::test::MakePoint(0.0, 0.0, 0.0);
  polygon.flux = 0.0;
  geometry.AddPolygon(polygon);

  mc3d::InnerBoundary boundary;
  boundary.SetGeometry(geometry);
  ASSERT_TRUE(boundary.AddPolygon(geometry,
                                  mc3d::test::MakePoint(-0.15, 0.0, 0.0), 0.5));
  ASSERT_FALSE(boundary.Empty());

  std::vector<Particle> particles{
      Particle(mc3d::test::MakePoint(-0.2, 0.0, 0.0),
               mc3d::test::MakePoint(1.0, 0.0, 0.0))};

  const Point initial_velocity = particles.front().velocity;
  constexpr double dt = 0.2;

  utils::SeedRandom(123u);
  EXPECT_EQ(boundary.BoundaryCondition(particles, dt), 0);

  const Point& final_position = particles.front().position;
  // Account for intentional offset (eps * 10.0 = 1e-4) added to prevent penetration
  EXPECT_NEAR(final_position.x, 0.0, 1e-3);
  EXPECT_NEAR(final_position.y, 0.0, 1e-9);
  EXPECT_NEAR(final_position.z, 0.0, 1e-9);

  mc3d::Polygon& collided_polygon = geometry.PolygonAt(0);
  const Point& normal = collided_polygon.GetNormal();
  const Point& final_velocity = particles.front().velocity;
  EXPECT_GT(normal * final_velocity, 0.0);

  const Point expected_force = Difference(initial_velocity, final_velocity);
  EXPECT_NEAR(collided_polygon.force.x, expected_force.x, kTolerance);
  EXPECT_NEAR(collided_polygon.force.y, expected_force.y, kTolerance);
  EXPECT_NEAR(collided_polygon.force.z, expected_force.z, kTolerance);
}

TEST(InnerBoundaryCollisionTest, ParticleCollidesWithBodyGeometry) {
  mc3d::test::GeometryAccessor geometry;
  geometry.CreateCube(-0.5, 1.0, 1.0, 1.0);

  mc3d::InnerBoundary boundary(geometry, mc3d::test::MakePoint(-0.6, 0.0, 0.0),
                               0.5);
  ASSERT_FALSE(boundary.Empty());

  std::vector<Particle> particles{
      Particle(mc3d::test::MakePoint(-0.6, 0.0, 0.0),
               mc3d::test::MakePoint(0.5, 0.2, 0.0))};

  const Point initial_velocity = particles.front().velocity;
  constexpr double dt = 0.2;

  utils::SeedRandom(456u);
  EXPECT_EQ(boundary.BoundaryCondition(particles, dt), 0);

  const Point& final_position = particles.front().position;
  // Account for intentional offset (eps * 10.0 = 1e-4) added to prevent penetration
  EXPECT_NEAR(final_position.x, -0.5, 1e-3);
  EXPECT_NEAR(final_position.y, initial_velocity.y * dt, 1e-9);
  EXPECT_NEAR(final_position.z, 0.0, 1e-9);

  const Point& final_velocity = particles.front().velocity;
  const auto polygons = geometry.RawPolygons();

  int polygons_with_force = 0;
  mc3d::Polygon* collided_polygon = nullptr;
  mc3d::Point accumulated_force(0.0, 0.0, 0.0);

  for (mc3d::Polygon* polygon : polygons) {
    if (polygon->force.Mod() > kTolerance) {
      polygons_with_force++;
      collided_polygon = polygon;
      accumulated_force += polygon->force;
    }
  }

  ASSERT_EQ(polygons_with_force, 1);
  ASSERT_NE(collided_polygon, nullptr);

  const Point expected_force = Difference(initial_velocity, final_velocity);
  EXPECT_NEAR(collided_polygon->force.x, expected_force.x, kTolerance);
  EXPECT_NEAR(collided_polygon->force.y, expected_force.y, kTolerance);
  EXPECT_NEAR(collided_polygon->force.z, expected_force.z, kTolerance);

  EXPECT_NEAR(accumulated_force.x, expected_force.x, kTolerance);
  EXPECT_NEAR(accumulated_force.y, expected_force.y, kTolerance);
  EXPECT_NEAR(accumulated_force.z, expected_force.z, kTolerance);

  const Point normal = collided_polygon->GetNormal();
  EXPECT_GT(normal * final_velocity, 0.0);
}
