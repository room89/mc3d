#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "geometry.h"
// Built-in cube and wedge are convex XY polygons extruded along Z.
// Clip in unit-cell coordinates, integrating disjoint outside pieces to avoid
// subtracting nearly equal volumes in thin gas layers. This query draws no RNG.
namespace mc3d {
namespace {
struct P {
  long double x, y;
};
using Poly = std::vector<P>;
long double Cross(P a, P b, P c) {
  return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}
Poly Clip(const Poly& p, P a, P b, bool inside) {
  Poly out;
  if (p.empty()) return out;
  for (size_t i = 0; i < p.size(); ++i) {
    P s = p[i], e = p[(i + 1) % p.size()];
    long double ds = Cross(a, b, s), de = Cross(a, b, e);
    if (!inside) {
      ds = -ds;
      de = -de;
    }
    if (ds >= 0) out.push_back(s);
    if ((ds < 0 && de > 0) || (ds > 0 && de < 0)) {
      long double t = ds / (ds - de);
      out.push_back({s.x + t * (e.x - s.x), s.y + t * (e.y - s.y)});
    }
  }
  return out;
}
struct Moment {
  long double area = 0, x = 0, y = 0;
};
Moment Integrate(const Poly& p) {
  Moment m;
  if (p.size() < 3) return m;
  P origin = p[0];
  for (size_t i = 1; i + 1 < p.size(); ++i) {
    long double a = Cross(origin, p[i], p[i + 1]) / 2;
    if (a < 0) throw std::runtime_error("Negative clipped polygon orientation");
    m.area += a;
    m.x += a * (origin.x + p[i].x + p[i + 1].x) / 3;
    m.y += a * (origin.y + p[i].y + p[i + 1].y) / 3;
  }
  return m;
}
}  // namespace
CellVolume Geometry::GasVolume(Point apex, Point size) const {
  if (!std::isfinite(apex.x) || !std::isfinite(apex.y) ||
      !std::isfinite(apex.z) || !(size.Volume() > 0) ||
      !(size.x > 0 && size.y > 0 && size.z > 0) ||
      !std::isfinite(size.Volume()) || apex.x + size.x == apex.x ||
      apex.y + size.y == apex.y || apex.z + size.z == apex.z)
    throw std::invalid_argument("Unresolvable cell dimensions");
  CellVolume result;
  result.full = size.Volume();
  if (PolygonCount() == 0) {
    result.gas = result.full;
    result.gas_center = apex + size / 2;
    return result;
  }
  if (!SupportsGasVolume())
    throw std::invalid_argument(
        "Gas volume supports only built-in cube and wedge; STL, cylinder and "
        "pyramid are unsupported");
  Poly vertices;
  long double zlo = INFINITY, zhi = -INFINITY;
  for (size_t i = 0; i < PolygonCount(); ++i) {
    auto& face = GetPolygon(i);
    for (Point v : {face.p1, face.p2, face.p3}) {
      if (!std::isfinite(v.x) || !std::isfinite(v.y) || !std::isfinite(v.z))
        throw std::invalid_argument("Nonfinite geometry coordinates");
      vertices.push_back({((long double)v.x - apex.x) / size.x,
                          ((long double)v.y - apex.y) / size.y});
      long double z = ((long double)v.z - apex.z) / size.z;
      zlo = std::min(zlo, z);
      zhi = std::max(zhi, z);
    }
  }
  std::sort(vertices.begin(), vertices.end(),
            [](P a, P b) { return a.x < b.x || (a.x == b.x && a.y < b.y); });
  vertices.erase(std::unique(vertices.begin(), vertices.end(),
                             [](P a, P b) { return a.x == b.x && a.y == b.y; }),
                 vertices.end());
  Poly hull;
  for (auto v : vertices) {
    while (hull.size() > 1 && Cross(hull[hull.size() - 2], hull.back(), v) <= 0)
      hull.pop_back();
    hull.push_back(v);
  }
  size_t lower = hull.size();
  for (auto it = vertices.rbegin() + 1; it != vertices.rend(); ++it) {
    while (hull.size() > lower &&
           Cross(hull[hull.size() - 2], hull.back(), *it) <= 0)
      hull.pop_back();
    hull.push_back(*it);
  }
  if (!hull.empty()) hull.pop_back();
  if (hull.size() < 3 || !(zhi > zlo))
    throw std::runtime_error("Degenerate supported body");
  long double a = std::clamp(zlo, 0.L, 1.L), b = std::clamp(zhi, 0.L, 1.L),
              depth = std::max(0.L, b - a);
  long double outside = a + (1 - b), gas = outside, gx = outside / 2,
              gy = outside / 2, gz = a * a / 2 + (1 - b) * (1 + b) / 2;
  Poly remaining{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  for (size_t i = 0; i < hull.size(); ++i) {
    P p = hull[i], q = hull[(i + 1) % hull.size()];
    auto m = Integrate(Clip(remaining, p, q, false));
    gas += m.area * depth;
    gx += m.x * depth;
    gy += m.y * depth;
    gz += m.area * depth * (a + b) / 2;
    remaining = Clip(remaining, p, q, true);
  }
  long double solid = Integrate(remaining).area * depth;
  if (!std::isfinite(gas) || !std::isfinite(solid) || gas < 0 || solid < 0 ||
      gas > 1 + 1e-12L || solid > 1 + 1e-12L ||
      std::abs(gas + solid - 1) > 1e-12L)
    throw std::runtime_error("Inconsistent clipped cell volumes");
  // Close the partition after the tolerance check; preserve every positive
  // gas fraction instead of flooring thin layers or subtracting from full.
  const long double partition = gas + solid;
  result.gas = double((gas / partition) * result.full);
  result.solid = double((solid / partition) * result.full);
  if (gas > 0 && result.gas == 0)
    throw std::runtime_error("Gas volume is not representable");
  if (result.gas > result.full || result.solid > result.full)
    throw std::runtime_error("Cell volume exceeds full volume");
  result.region = gas == 0     ? CellRegion::Solid
                  : solid == 0 ? CellRegion::Gas
                               : CellRegion::Cut;
  if (gas > 0)
    result.gas_center = Point(apex.x + double(gx / gas) * size.x,
                              apex.y + double(gy / gas) * size.y,
                              apex.z + double(gz / gas) * size.z);
  return result;
}
}  // namespace mc3d
