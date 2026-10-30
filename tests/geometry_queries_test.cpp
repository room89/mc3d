// Frozen query algorithm from bfb108f: oracle independent of prepared kernels.
#include "geometry.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <future>
#include <limits>
#include <random>

namespace {
using namespace mc3d;
// Returns the fraction of a ray needed to reach a triangle. The ray direction
// can be a velocity step or an arbitrary direction for point-in-mesh queries.
std::optional<double> ReferenceRayTriangle(Point start, Point direction,
                                  const Polygon& polygon) {
  const Point e1 = polygon.GetP2() - polygon.GetP1();
  const Point e2 = polygon.GetP3() - polygon.GetP1();
  const Point p = direction.Cross(e2);
  const double det = e1 * p;
  if (std::abs(det) <= 1e-14 * e1.Mod() * e2.Mod() * direction.Mod()) {
    return std::nullopt;
  }
  const double inverse_det = 1.0 / det;
  Point from_vertex = start - polygon.GetP1();
  const double u = (from_vertex * p) * inverse_det;
  if (u < -1e-10 || u > 1.0 + 1e-10) return std::nullopt;
  const Point q = from_vertex.Cross(e1);
  const double v = (direction * q) * inverse_det;
  if (v < -1e-10 || u + v > 1.0 + 1e-10) return std::nullopt;
  return (e2 * q) * inverse_det;
}


class ReferenceGeometry : public Geometry {
 public:
  bool ReferenceIsInnerPoint(Point point) const;
  std::optional<SurfaceHit> ReferenceFirstIntersection(Point start, Point displacement);
  std::optional<Point> ReferenceExteriorPoint(Point point) const;
  void Triangle(Point a, Point b, Point c, Point normal) {
    InvalidateQueries();
    poligons.push_back(std::make_unique<Polygon>(a,b,c,normal));
  }
  void AddTriangle(Point a, Point b, Point c) {
    InvalidateQueries();
    poligons.push_back(std::make_unique<Polygon>(a,b,c));
  }
};
bool ReferenceGeometry::ReferenceIsInnerPoint(Point test_point) const {
  // A single surface triangle has no enclosed volume; at least four faces
  // are needed for a closed three-dimensional body.
  if (poligons.size() < 4) {
    return false;
  }
  const auto [lower, upper] = Bounds();
  // Expand beyond the ray/triangle barycentric tolerance to retain near-surface behavior.
  const double margin = std::max(SurfaceTolerance(), 1e-9 * (upper - lower).Mod());
  if (test_point.x < lower.x - margin || test_point.x > upper.x + margin ||
      test_point.y < lower.y - margin || test_point.y > upper.y + margin ||
      test_point.z < lower.z - margin || test_point.z > upper.z + margin) return false;
  // A parity test works for closed concave meshes as well as convex bodies.
  // Three non-axis-aligned rays avoid most edge and vertex degeneracies.
  const std::array<Point, 3> directions{
      {Point(1, 0.371, 0.529), Point(0.271, 1, 0.419), Point(0.379, 0.233, 1)}};
  const double tolerance =
      std::max(SurfaceTolerance(),
               16 * std::numeric_limits<double>::epsilon() *
                   std::max({std::abs(test_point.x), std::abs(test_point.y),
                             std::abs(test_point.z)}));
  int inside_votes = 0;
  // Queries do not recurse; each worker can safely reuse its own scratch storage.
  thread_local std::vector<double> distances;
  distances.reserve(poligons.size());
  for (const Point& direction : directions) {
    distances.clear();
    for (const auto& polygon : poligons) {
      const auto fraction = ReferenceRayTriangle(test_point, direction, *polygon);
      if (!fraction) continue;
      if (std::abs(*fraction) * direction.Mod() <= tolerance) return false;
      if (*fraction > 0) distances.push_back(*fraction);
    }
    std::sort(distances.begin(), distances.end());
    std::size_t crossings = 0;
    double previous = -std::numeric_limits<double>::infinity();
    for (double distance : distances) {
      if ((distance - previous) * direction.Mod() > tolerance) {
        ++crossings;
        previous = distance;
      }
    }
    inside_votes += crossings % 2;
  }
  return inside_votes >= 2;
}

std::optional<ReferenceGeometry::SurfaceHit> ReferenceGeometry::ReferenceFirstIntersection(
    Point start, Point displacement) {
  std::optional<SurfaceHit> first;
  for (const auto& polygon : poligons) {
    const auto fraction = ReferenceRayTriangle(start, displacement, *polygon);
    // A positive fraction is a real hit even for a very long segment.
    if (!fraction || *fraction < 0 || *fraction > 1.0) continue;
    if (*fraction == 0) {
      const double distance = displacement.Mod();
      if (distance == 0 ||
          !ReferenceIsInnerPoint(start +
                        displacement *
                            std::min(1.0, 8 * SurfaceTolerance() / distance))) {
        continue;
      }
    }
    if (!first || *fraction < first->fraction) {
      first = SurfaceHit{polygon.get(), *fraction};
    }
  }
  return first;
}

std::optional<Point> ReferenceGeometry::ReferenceExteriorPoint(Point interior) const {
  if (!ReferenceIsInnerPoint(interior)) return interior;
  const std::array<Point, 6> directions{{Point(1, 0, 0), Point(-1, 0, 0),
                                         Point(0, 1, 0), Point(0, -1, 0),
                                         Point(0, 0, 1), Point(0, 0, -1)}};
  double nearest = std::numeric_limits<double>::infinity();
  Point exit;
  Point direction_to_exit;
  for (const Point& direction : directions) {
    for (const auto& polygon : poligons) {
      const auto distance = ReferenceRayTriangle(interior, direction, *polygon);
      if (distance && *distance > 0 && *distance < nearest) {
        nearest = *distance;
        exit = interior + direction * *distance;
        direction_to_exit = direction;
      }
    }
  }
  if (!std::isfinite(nearest)) return std::nullopt;
  // Use a geometry-scaled offset that is representable at these coordinates;
  // check the result before accepting the recovery.
  double offset = 8 * SurfaceTolerance();
  for (int attempt = 0; attempt < 12; ++attempt) {
    Point candidate = exit + direction_to_exit * offset;
    if (!ReferenceIsInnerPoint(candidate)) return candidate;
    offset *= 2;
  }
  return std::nullopt;
}


} // namespace

namespace {
std::unique_ptr<ReferenceGeometry> MakeConcavePrism() {
  auto geometry = std::make_unique<ReferenceGeometry>();
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

void Compare(ReferenceGeometry& body, Point start, Point displacement) {
  EXPECT_EQ(body.IsInnerPoint(start), body.ReferenceIsInnerPoint(start));
  const auto actual = body.FirstIntersection(start, displacement);
  const auto expected = body.ReferenceFirstIntersection(start, displacement);
  ASSERT_EQ(actual.has_value(), expected.has_value());
  if (actual) {
    EXPECT_EQ(actual->polygon, expected->polygon);
    EXPECT_EQ(actual->fraction, expected->fraction);
  }
  const auto exit = body.ExteriorPoint(start);
  const auto reference_exit = body.ReferenceExteriorPoint(start);
  ASSERT_EQ(exit.has_value(), reference_exit.has_value());
  if (exit) {
    EXPECT_EQ(exit->x, reference_exit->x);
    EXPECT_EQ(exit->y, reference_exit->y);
    EXPECT_EQ(exit->z, reference_exit->z);
  }
}

void CheckQueries(ReferenceGeometry& body) {
  const auto [lo, hi] = body.Bounds();
  const Point size = hi - lo;
  const Point center = (lo + hi) * .5;
  for (Point start : {lo, hi, center, Point(hi.x, center.y, center.z)})
    for (Point direction : {Point(0,0,0), size, size * -1., Point(1e12,0,0), Point(1,1e-15,0)})
      Compare(body, start, direction);
  for (std::size_t i = 0; i < body.PolygonCount(); ++i) {
    const auto& triangle = body.GetPolygon(i);
    const auto a = triangle.GetP1(), b = triangle.GetP2(), c = triangle.GetP3();
    Point normal = (b-a).Cross(c-a);
    if (normal.Mod() > 0) normal /= normal.Mod();
    for (Point point : {a, (a+b)*.5, (a+b+c)/3.})
      for (double offset : {-2., -.5, 0., .5, 2.})
        Compare(body, point + normal * (offset * body.SurfaceTolerance()), normal);
  }
  std::mt19937 rng(9841);
  std::uniform_real_distribution<double> u(-.5,1.5);
  for (int i = 0; i < 150; ++i) {
    const Point start = lo + Point(size.x*u(rng),size.y*u(rng),size.z*u(rng));
    const Point end = lo + Point(size.x*u(rng),size.y*u(rng),size.z*u(rng));
    Compare(body, start, end-start);
  }
}
}

TEST(PreparedGeometryTest, MatchesFrozenOracleAcrossMeshesAndScales) {
  for (int shape = 0; shape < 4; ++shape) {
    auto body = shape == 3 ? MakeConcavePrism() : std::make_unique<ReferenceGeometry>();
    if (shape == 0) body->CreateCube(-.5,1,1,1);
    if (shape == 1) body->CreateWedge();
    if (shape == 2) body->CreateCylinder(-.5,.5,1,32);
    for (int transform = 0; transform < 3; ++transform) {
      SCOPED_TRACE(shape);
      SCOPED_TRACE(transform);
      if (transform == 1) body->Scale(1e-6);
      if (transform == 2) { body->Scale(1e6); body->Move({1e8,-1e8,1e8}); }
      CheckQueries(*body);
      body->PrepareQueries();
      CheckQueries(*body);
      body->PrepareQueries();
      CheckQueries(*body);
    }
  }
}

TEST(PreparedGeometryTest, InvalidatesAfterEverySupportedEdit) {
  ReferenceGeometry body;
  body.CreateCube(-.5,1,1,1);
  for (int edit = 0; edit < 7; ++edit) {
    body.PrepareQueries();
    switch (edit) {
      case 0: body.Move({2,-3,1}); break;
      case 1: body.Scale(.5); break;
      case 2: body.Fragment(.4); break;
      case 3: body.FixPolygons(); break;
      case 4: body.ReverseNormals(); break;
      case 5: body.CreateWedge(); break;
      case 6: body.AddTriangle({5,0,0},{5,1,0},{5,0,1}); break;
    }
    CheckQueries(body);
    body.PrepareQueries();
    CheckQueries(body);
  }
}

TEST(PreparedGeometryTest, EmptyDegenerateAndCoincidentTriangles) {
  ReferenceGeometry body;
  body.PrepareQueries();
  Compare(body,{0,0,0},{1,0,0});
  body.AddTriangle({0,0,0},{0,0,0},{0,0,0});
  body.AddTriangle({0,0,0},{0,1,0},{0,0,1});
  body.AddTriangle({0,0,0},{0,1,0},{0,0,1});
  body.AddTriangle({0,0,0},{0,1,0},{0,0,1});
  for (bool prepared : {false,true}) {
    if (prepared) body.PrepareQueries();
    CheckQueries(body);
    auto hit = body.FirstIntersection({-1,.2,.2},{2,0,0});
    ASSERT_TRUE(hit);
    EXPECT_EQ(hit->polygon, &body.GetPolygon(1));
  }
}

TEST(PreparedGeometryTest, ConcurrentReadOnlyQueriesMatchOracle) {
  auto body = MakeConcavePrism();
  body->PrepareQueries();
  std::vector<std::future<void>> workers;
  for (int i = 0; i < 4; ++i)
    workers.push_back(std::async(std::launch::async, [&] { CheckQueries(*body); }));
  for (auto& worker : workers) worker.get();
}
