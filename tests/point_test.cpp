#include <gtest/gtest.h>
#include <point.h>

#include <utils/utils.hpp>

namespace {

constexpr double kTolerance = 1e-12;

}  // namespace

TEST(PointTest, NormalizeProducesUnitVector) {
  mc3d::Point point(3.0, 4.0, 0.0);
  mc3d::Point& result = point.Normalize();

  EXPECT_EQ(&result, &point);
  EXPECT_NEAR(point.x, 0.6, kTolerance);
  EXPECT_NEAR(point.y, 0.8, kTolerance);
  EXPECT_NEAR(point.z, 0.0, kTolerance);
  EXPECT_NEAR(point.Mod(), 1.0, kTolerance);
}

TEST(PointTest, CrossProductFollowsRightHandRule) {
  mc3d::Point a(1.0, 0.0, 0.0);
  mc3d::Point b(0.0, 1.0, 0.0);

  mc3d::Point cross_ab = a.Cross(b);
  EXPECT_NEAR(cross_ab.x, 0.0, kTolerance);
  EXPECT_NEAR(cross_ab.y, 0.0, kTolerance);
  EXPECT_NEAR(cross_ab.z, 1.0, kTolerance);

  mc3d::Point cross_ba = b.Cross(a);
  EXPECT_NEAR(cross_ba.x, 0.0, kTolerance);
  EXPECT_NEAR(cross_ba.y, 0.0, kTolerance);
  EXPECT_NEAR(cross_ba.z, -1.0, kTolerance);
}

TEST(PointTest, VolumeMultipliesComponents) {
  const mc3d::Point box(2.5, 4.0, 1.5);
  EXPECT_DOUBLE_EQ(box.Volume(), 15.0);
}

TEST(PointRandomTest, StaticRandomPointRespectsBounds) {
  const mc3d::Point bounds(2.0, 3.0, 4.0);
  const mc3d::Point sample = mc3d::Point::RandomPoint(bounds);

  EXPECT_GE(sample.x, 0.0);
  EXPECT_LE(sample.x, bounds.x);
  EXPECT_GE(sample.y, 0.0);
  EXPECT_LE(sample.y, bounds.y);
  EXPECT_GE(sample.z, 0.0);
  EXPECT_LE(sample.z, bounds.z);
}

TEST(PointRandomTest, MemberRandomPointRespectsBounds) {
  mc3d::Point bounds(5.0, 0.0, 1.0);
  const mc3d::Point sample = bounds.RandomPoint();

  EXPECT_GE(sample.x, 0.0);
  EXPECT_LE(sample.x, bounds.x);
  EXPECT_DOUBLE_EQ(sample.y, 0.0);
  EXPECT_GE(sample.z, 0.0);
  EXPECT_LE(sample.z, bounds.z);
}
