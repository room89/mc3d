#include "solver_setup.h"

#include <gtest/gtest.h>
#include <stdexcept>

#include "free_boundary.h"
#include "giper_free_boundary.h"
#include "mirror_boundary.h"
#include "pereodic_boundary.h"
#include "poligon.h"

namespace {

TEST(SimulationConfigTest, ParsesAndValidatesWallTemperature) {
  char program[] = "solver";
  char option[] = "--wall-temperature";
  char valid[] = "10";
  char invalid[] = "0";
  char* args[] = {program, option, valid};
  EXPECT_DOUBLE_EQ(mc3d::LoadSimulationConfig(1, args).wall_temperature, 1);
  EXPECT_DOUBLE_EQ(mc3d::LoadSimulationConfig(3, args).wall_temperature, 10);
  args[2] = invalid;
  EXPECT_THROW(mc3d::LoadSimulationConfig(3, args), std::invalid_argument);
}

constexpr double kTolerance = 1e-9;

struct DescriptorExpectation {
  mc3d::BoundaryFace face;
  mc3d::Point expected_position;
  mc3d::Point expected_normal;
  mc3d::Point expected_mixing;
};

class GeometryBuilder : public mc3d::Geometry {
 public:
  void AddPolygon(const mc3d::Point& p1, const mc3d::Point& p2,
                  const mc3d::Point& p3, const mc3d::Point& normal) {
    poligons.push_back(std::make_unique<mc3d::Polygon>(p1, p2, p3, normal));
  }

  mc3d::Polygon& PolygonAt(std::size_t index) { return *poligons.at(index); }
};

void ExpectPointEqual(const mc3d::Point& actual, const mc3d::Point& expected) {
  EXPECT_NEAR(actual.x, expected.x, kTolerance);
  EXPECT_NEAR(actual.y, expected.y, kTolerance);
  EXPECT_NEAR(actual.z, expected.z, kTolerance);
}

double SignedOrientation(const mc3d::Polygon& polygon) {
  mc3d::Point edge1 = polygon.GetP2() - polygon.GetP1();
  mc3d::Point edge2 = polygon.GetP3() - polygon.GetP2();
  return edge1.Cross(edge2) * polygon.GetNormal();
}

}  // namespace

TEST(SolverSetupTest, MakeBoundaryDescriptorProducesExpectedValues) {
  mc3d::SimulationConfig cfg;
  cfg.Lx = 2;
  cfg.Ly = 3;
  cfg.Lz = 5;

  const DescriptorExpectation expectations[] = {
      {mc3d::BoundaryFace::XNeg, mc3d::Point(-0.5 * cfg.Lx, 0.0, 0.0),
       mc3d::Point(-1.0, 0.0, 0.0), mc3d::Point(cfg.Lx, 0.0, 0.0)},
      {mc3d::BoundaryFace::XPos, mc3d::Point(0.5 * cfg.Lx, 0.0, 0.0),
       mc3d::Point(1.0, 0.0, 0.0), mc3d::Point(-cfg.Lx, 0.0, 0.0)},
      {mc3d::BoundaryFace::YNeg, mc3d::Point(0.0, -0.5 * cfg.Ly, 0.0),
       mc3d::Point(0.0, -1.0, 0.0), mc3d::Point(0.0, cfg.Ly, 0.0)},
      {mc3d::BoundaryFace::YPos, mc3d::Point(0.0, 0.5 * cfg.Ly, 0.0),
       mc3d::Point(0.0, 1.0, 0.0), mc3d::Point(0.0, -cfg.Ly, 0.0)},
      {mc3d::BoundaryFace::ZNeg, mc3d::Point(0.0, 0.0, -0.5 * cfg.Lz),
       mc3d::Point(0.0, 0.0, -1.0), mc3d::Point(0.0, 0.0, cfg.Lz)},
      {mc3d::BoundaryFace::ZPos, mc3d::Point(0.0, 0.0, 0.5 * cfg.Lz),
       mc3d::Point(0.0, 0.0, 1.0), mc3d::Point(0.0, 0.0, -cfg.Lz)}};

  for (const auto& expectation : expectations) {
    const auto descriptor = mc3d::MakeBoundaryDescriptor(expectation.face, cfg);
    ExpectPointEqual(descriptor.position, expectation.expected_position);
    ExpectPointEqual(descriptor.normal, expectation.expected_normal);
    ExpectPointEqual(descriptor.mixing, expectation.expected_mixing);
  }
}

TEST(SolverSetupTest, MakeBoundaryReturnsCorrectSubtype) {
  mc3d::SimulationConfig cfg;
  const auto descriptor =
      mc3d::MakeBoundaryDescriptor(mc3d::BoundaryFace::XNeg, cfg);

  {
    auto boundary =
        mc3d::MakeBoundary(descriptor, mc3d::BoundaryType::Mirror, cfg);
    ASSERT_NE(boundary, nullptr);
    auto* typed = dynamic_cast<mc3d::MirrorBoundary*>(boundary.get());
    ASSERT_NE(typed, nullptr);
    ExpectPointEqual(typed->GetNormal(), mc3d::Point(-1.0, 0.0, 0.0));
    ExpectPointEqual(typed->GetPosition(),
                     mc3d::Point(-0.5 * cfg.Lx, 0.0, 0.0));
  }

  {
    auto boundary =
        mc3d::MakeBoundary(descriptor, mc3d::BoundaryType::Periodic, cfg);
    ASSERT_NE(boundary, nullptr);
    auto* typed = dynamic_cast<mc3d::PeriodicBoundary*>(boundary.get());
    ASSERT_NE(typed, nullptr);
    ExpectPointEqual(typed->GetMixing(), mc3d::Point(cfg.Lx, 0.0, 0.0));
  }

  {
    auto boundary =
        mc3d::MakeBoundary(descriptor, mc3d::BoundaryType::Free, cfg);
    ASSERT_NE(boundary, nullptr);
    EXPECT_NE(dynamic_cast<mc3d::FreeBoundary*>(boundary.get()), nullptr);
  }

  {
    auto boundary =
        mc3d::MakeBoundary(descriptor, mc3d::BoundaryType::HyperFree, cfg);
    ASSERT_NE(boundary, nullptr);
    EXPECT_NE(dynamic_cast<mc3d::HyperFreeBoundary*>(boundary.get()), nullptr);
  }

  auto none = mc3d::MakeBoundary(descriptor, mc3d::BoundaryType::None, cfg);
  EXPECT_EQ(none, nullptr);
}

TEST(SolverSetupTest, BuildGeometryCreatesWedgeWhenRequested) {
  mc3d::SimulationConfig cfg;
  cfg.geometry_type = "wedge";

  auto geometry = mc3d::BuildGeometry(cfg);
  ASSERT_NE(geometry, nullptr);
  EXPECT_EQ(geometry->PolygonCount(), 8U);
}

TEST(SolverSetupTest, BuildGeometryCreatesCubeWithConfiguredDimensions) {
  mc3d::SimulationConfig cfg;
  cfg.geometry_type = "cube";
  cfg.geometry_cube_x = -0.3;
  cfg.geometry_cube_length = 0.6;
  cfg.geometry_cube_width = 0.4;
  cfg.geometry_cube_height = 0.8;

  auto geometry = mc3d::BuildGeometry(cfg);
  ASSERT_NE(geometry, nullptr);
  EXPECT_EQ(geometry->PolygonCount(), 12U);

  const auto size = geometry->Size();
  ExpectPointEqual(size.first, mc3d::Point(cfg.geometry_cube_x,
                                           -cfg.geometry_cube_height / 2,
                                           -cfg.geometry_cube_width / 2));
  ExpectPointEqual(
      size.second,
      mc3d::Point(cfg.geometry_cube_x + cfg.geometry_cube_length,
                  cfg.geometry_cube_height / 2, cfg.geometry_cube_width / 2));

  const auto mass_center = geometry->MassCenter();
  ExpectPointEqual(mass_center, mc3d::Point(cfg.geometry_cube_x +
                                                cfg.geometry_cube_length / 2,
                                            0.0, 0.0));
}

TEST(SolverSetupTest, BuildGeometryAppliesPostProcessingOptions) {
  mc3d::SimulationConfig cfg;
  cfg.geometry_type = "cube";
  cfg.geometry_cube_x = -1.0;
  cfg.geometry_cube_length = 2.0;
  cfg.geometry_cube_width = 0.4;
  cfg.geometry_cube_height = 0.8;
  cfg.geometry_reverse_normals = true;
  cfg.geometry_scale = 2.0;
  cfg.geometry_move_x = 0.5;
  cfg.geometry_move_y = -0.25;
  cfg.geometry_move_z = 0.1;

  auto geometry = mc3d::BuildGeometry(cfg);
  ASSERT_NE(geometry, nullptr);

  const auto size = geometry->Size();
  mc3d::Point min_expected(cfg.geometry_cube_x, -cfg.geometry_cube_height / 2,
                           -cfg.geometry_cube_width / 2);
  mc3d::Point max_expected(cfg.geometry_cube_x + cfg.geometry_cube_length,
                           cfg.geometry_cube_height / 2,
                           cfg.geometry_cube_width / 2);
  min_expected *= *cfg.geometry_scale;
  max_expected *= *cfg.geometry_scale;
  const mc3d::Point delta(cfg.geometry_move_x.value(),
                          cfg.geometry_move_y.value(),
                          cfg.geometry_move_z.value());
  min_expected += delta;
  max_expected += delta;

  ExpectPointEqual(size.first, min_expected);
  ExpectPointEqual(size.second, max_expected);

  const auto& polygon = geometry->GetPolygon(0);
  ExpectPointEqual(polygon.GetNormal(), mc3d::Point(1.0, 0.0, 0.0));
}

TEST(SolverSetupTest, ApplyGeometryPostProcessingFixesPolygonOrientation) {
  mc3d::SimulationConfig cfg;
  cfg.geometry_fix_polygons = true;

  GeometryBuilder geometry;
  geometry.AddPolygon(mc3d::Point(0, 0, 0), mc3d::Point(1, 0, 0),
                      mc3d::Point(0, 1, 0), mc3d::Point(0, 0, -1));

  mc3d::ApplyGeometryPostProcessing(cfg, geometry);

  ASSERT_EQ(geometry.PolygonCount(), 1U);
  const auto& polygon = geometry.PolygonAt(0);
  EXPECT_GT(SignedOrientation(polygon), 0.0);
}

TEST(SolverSetupTest, ApplyGeometryPostProcessingFragmentsGeometry) {
  mc3d::SimulationConfig cfg;
  cfg.geometry_fragment_length = 0.25;

  GeometryBuilder geometry;
  geometry.AddPolygon(mc3d::Point(0, 0, 0), mc3d::Point(1, 0, 0),
                      mc3d::Point(0, 1, 0), mc3d::Point(0, 0, 1));
  const auto before_count = geometry.PolygonCount();

  mc3d::ApplyGeometryPostProcessing(cfg, geometry);

  EXPECT_GT(geometry.PolygonCount(), before_count);
}
