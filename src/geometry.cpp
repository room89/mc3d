#include "geometry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace mc3d {
namespace {

// Returns the fraction of a ray needed to reach a triangle. The ray direction
// can be a velocity step or an arbitrary direction for point-in-mesh queries.
std::optional<double> RayTriangle(Point start, Point direction,
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

}  // namespace
Geometry::Geometry() {}

Geometry::Geometry(const char* file_name) {
  std::cout << endl
            << "...reading Geometry from file " << file_name << "... " << endl;
  string word;
  vector<string> v;
  std::ifstream file(file_name);
  if (!file.is_open()) {
    std::cout << endl << "Can't open file!!!" << endl;
    throw unusual_situations::Exception(
        mc3d::unusual_situations::exit_code::BOUNDARY_CONSTRUCTOR_ERROR);
  }
  while (file >> word) {
    v.push_back(word);
  }
  if (v[0] == "solid" && v[v.size() - 1] == "endsolid") {
    unsigned int i = 1;
    unsigned int n = 0;
    if (v[1] != "facet") {
      std::cout << endl << "Solid: " << v[1] << endl;
      i++;
    }
    while (i < v.size()) {
      if (v[i] == "facet" && v[i + 1] == "normal" && v[i + 6] == "loop" &&
          v[i + 19] == "endloop" && v[i + 20] == "endfacet") {
        if (v[i] == "endsolid") break;
        auto temp = std::make_unique<Polygon>();
        mc3d::Point p1, p2, p3, normal;

        p1.x = std::atof(v[i + 8].c_str());
        p1.y = std::atof(v[i + 9].c_str());
        p1.z = std::atof(v[i + 10].c_str());

        p2.x = std::atof(v[i + 12].c_str());
        p2.y = std::atof(v[i + 13].c_str());
        p2.z = std::atof(v[i + 14].c_str());

        p3.x = std::atof(v[i + 16].c_str());
        p3.y = std::atof(v[i + 17].c_str());
        p3.z = std::atof(v[i + 18].c_str());

        normal.x = std::atof(v[i + 2].c_str());
        normal.y = std::atof(v[i + 3].c_str());
        normal.z = std::atof(v[i + 4].c_str());

        normal.Normalize();

        temp->Set(p1, p2, p3, normal);

        poligons.push_back(std::move(temp));

        i += 21;
      } else {
        i++;
      }
    }
  }

  file.close();
}

Geometry::~Geometry() = default;

void Geometry::WriteGeometryFile(const char* file_name, const char* body_name) {
  std::ofstream file(file_name);

  if (!file.is_open()) return;

  file << "solid " << body_name << std::endl;

  for (const auto& polygon : poligons) {
    file << "facet normal " << std::setw(15) << polygon->GetNormal().x
         << std::setw(15) << polygon->GetNormal().y << std::setw(15)
         << polygon->GetNormal().z << std::endl;
    file << "\touter loop" << std::endl;
    file << "\t\tvertex" << std::setw(15) << polygon->GetP1().x << std::setw(15)
         << polygon->GetP1().y << std::setw(15) << polygon->GetP1().z
         << std::endl;
    file << "\t\tvertex" << std::setw(15) << polygon->GetP2().x << std::setw(15)
         << polygon->GetP2().y << std::setw(15) << polygon->GetP2().z
         << std::endl;
    file << "\t\tvertex" << std::setw(15) << polygon->GetP3().x << std::setw(15)
         << polygon->GetP3().y << std::setw(15) << polygon->GetP3().z
         << std::endl;
    file << "\tendloop" << std::endl << "endfacet" << std::endl;
  }

  file << "endsolid" << body_name << std::endl;

  file.close();
}

Point Geometry::MassCenter() {
  Point mc(0, 0, 0);

  double Stot = 0;

  for (const auto& polygon : poligons) {
    Point a = polygon->GetP2() - polygon->GetP1();
    Point b = polygon->GetP3() - polygon->GetP2();

    double S = 0.5 * a.Cross(b).Mod();

    Stot += S;
    mc += S * (polygon->GetP1() + polygon->GetP2() + polygon->GetP3()) / 3;
  }

  if (Stot > 0) {
    mc /= Stot;
  }

  return mc;
}

void Geometry::Move(Point a) {
  for (auto& polygon : poligons) {
    polygon->Move(a);
  }
}

void Geometry::Scale(double e) {
  for (auto& polygon : poligons) {
    polygon->Scale(e);
  }
}

std::size_t Geometry::PolygonCount() const { return poligons.size(); }

const Polygon& Geometry::GetPolygon(std::size_t index) const {
  return *poligons.at(index);
}

bool Geometry::IsInnerPoint(Point test_point) const {
  // A single surface triangle has no enclosed volume; at least four faces
  // are needed for a closed three-dimensional body.
  if (poligons.size() < 4) {
    return false;
  }
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
  for (const Point& direction : directions) {
    std::vector<double> distances;
    distances.reserve(poligons.size());
    for (const auto& polygon : poligons) {
      const auto fraction = RayTriangle(test_point, direction, *polygon);
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

std::optional<Geometry::SurfaceHit> Geometry::FirstIntersection(
    Point start, Point displacement) {
  std::optional<SurfaceHit> first;
  for (const auto& polygon : poligons) {
    const auto fraction = RayTriangle(start, displacement, *polygon);
    // A positive fraction is a real hit even for a very long segment.
    if (!fraction || *fraction < 0 || *fraction > 1.0) continue;
    if (*fraction == 0) {
      const double distance = displacement.Mod();
      if (distance == 0 ||
          !IsInnerPoint(start +
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

std::optional<Point> Geometry::ExteriorPoint(Point interior) const {
  if (!IsInnerPoint(interior)) return interior;
  const std::array<Point, 6> directions{{Point(1, 0, 0), Point(-1, 0, 0),
                                         Point(0, 1, 0), Point(0, -1, 0),
                                         Point(0, 0, 1), Point(0, 0, -1)}};
  double nearest = std::numeric_limits<double>::infinity();
  Point exit;
  Point direction_to_exit;
  for (const Point& direction : directions) {
    for (const auto& polygon : poligons) {
      const auto distance = RayTriangle(interior, direction, *polygon);
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
    if (!IsInnerPoint(candidate)) return candidate;
    offset *= 2;
  }
  return std::nullopt;
}

void Geometry::AccumulateForce(Polygon& polygon, Point impulse) {
  std::lock_guard<std::mutex> lock(force_mutex_);
  polygon.force += impulse;
}

std::pair<Point, Point> Geometry::Bounds() const {
  Point lower(std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity());
  Point upper(-std::numeric_limits<double>::infinity(),
              -std::numeric_limits<double>::infinity(),
              -std::numeric_limits<double>::infinity());
  for (const auto& polygon : poligons) {
    for (const Point& vertex :
         {polygon->GetP1(), polygon->GetP2(), polygon->GetP3()}) {
      lower.x = std::min(lower.x, vertex.x);
      lower.y = std::min(lower.y, vertex.y);
      lower.z = std::min(lower.z, vertex.z);
      upper.x = std::max(upper.x, vertex.x);
      upper.y = std::max(upper.y, vertex.y);
      upper.z = std::max(upper.z, vertex.z);
    }
  }
  return {lower, upper};
}

double Geometry::SurfaceTolerance() const {
  if (poligons.empty()) return 0;
  const auto [lower, upper] = Bounds();
  const Point extent = upper - lower;
  const double scale = std::max({extent.x, extent.y, extent.z});
  const double coordinate =
      std::max({std::abs(lower.x), std::abs(lower.y), std::abs(lower.z),
                std::abs(upper.x), std::abs(upper.y), std::abs(upper.z)});
  return std::max(1e-12 * scale,
                  16 * std::numeric_limits<double>::epsilon() * coordinate);
}

void Geometry::CreateWedge(double x, double width, double length,
                           double alpha) {
  Point vertices[6];

  vertices[0].Set(x, 0, -(width / 2));
  vertices[1].Set(x + length, length * sin(alpha), -width / 2);
  vertices[2].Set(x + length, -length * sin(alpha), -width / 2);

  vertices[3].Set(x, 0, width / 2);
  vertices[4].Set(x + length, length * sin(alpha), width / 2);
  vertices[5].Set(x + length, -length * sin(alpha), width / 2);

  auto add_polygon = [&](const Point& p1, const Point& p2, const Point& p3,
                         const Point& normal) {
    auto polygon = std::make_unique<Polygon>(p1, p2, p3, normal);
    polygon->flux = 0;
    polygon->force = Point(0, 0, 0);
    poligons.push_back(std::move(polygon));
  };

  add_polygon(vertices[0], vertices[1], vertices[2],
              -1 * (vertices[0] - vertices[1])
                       .Cross(vertices[2] - vertices[1])
                       .Normalize());
  add_polygon(vertices[3], vertices[4], vertices[5],
              -1 * (vertices[4] - vertices[3])
                       .Cross(vertices[5] - vertices[4])
                       .Normalize());
  add_polygon(vertices[0], vertices[1], vertices[3],
              -1 * (vertices[1] - vertices[0])
                       .Cross(vertices[3] - vertices[0])
                       .Normalize());
  add_polygon(vertices[1], vertices[4], vertices[3],
              -1 * (vertices[4] - vertices[1])
                       .Cross(vertices[3] - vertices[1])
                       .Normalize());
  add_polygon(vertices[0], vertices[2], vertices[3],
              -1 * (vertices[3] - vertices[0])
                       .Cross(vertices[2] - vertices[0])
                       .Normalize());
  add_polygon(vertices[2], vertices[3], vertices[5],
              -1 * (vertices[3] - vertices[2])
                       .Cross(vertices[5] - vertices[2])
                       .Normalize());
  add_polygon(vertices[1], vertices[2], vertices[5],
              -1 * (vertices[2] - vertices[1])
                       .Cross(vertices[4] - vertices[1])
                       .Normalize());
  add_polygon(vertices[1], vertices[5], vertices[4],
              -1 * (vertices[5] - vertices[1])
                       .Cross(vertices[4] - vertices[1])
                       .Normalize());
}

void Geometry::Fragment(double Lmax) {
  double max_Lmax = 10000.;
  while (max_Lmax > Lmax) {
    max_Lmax = 0;
    auto poligon_iter = poligons.begin();

    while (poligon_iter != poligons.end()) {
      if ((*poligon_iter)->GetLmax() > Lmax) {
        max_Lmax > (*poligon_iter)->GetLmax()
            ? max_Lmax = max_Lmax
            : max_Lmax = (*poligon_iter)->GetLmax();
        auto new_poligon = (*poligon_iter)->Divide();
        poligon_iter = poligons.erase(poligon_iter);
        poligons.push_back(std::move(new_poligon.first));
        poligons.push_back(std::move(new_poligon.second));
        poligon_iter = poligons.begin();
        continue;
      }

      poligon_iter++;
    }
  }
}

std::pair<Point, Point> Geometry::Size() {
  Point p1(0, 0, 0), p2(0, 0, 0);

  auto poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    p1.x = std::max(p1.x, (*poligon_iter)->p1.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p1.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p1.z);

    p1.x = std::max(p1.x, (*poligon_iter)->p2.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p2.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p2.z);

    p1.x = std::max(p1.x, (*poligon_iter)->p3.x);
    p1.y = std::max(p1.y, (*poligon_iter)->p3.y);
    p1.z = std::max(p1.z, (*poligon_iter)->p3.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p1.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p1.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p1.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p2.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p2.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p2.z);

    p2.x = std::min(p2.x, (*poligon_iter)->p3.x);
    p2.y = std::min(p2.y, (*poligon_iter)->p3.y);
    p2.z = std::min(p2.z, (*poligon_iter)->p3.z);

    poligon_iter++;
  }

  return pair<Point, Point>(p2, p1);
}

void Geometry::CreatePyramid(double x, double width, double length, double H) {}

double Geometry::DistanceToPoint(Point p) {
  double dist = 100000.;

  auto poligon_iter = poligons.begin();

  while (poligon_iter != poligons.end()) {
    dist = std::min(dist, (*poligon_iter)->DistanceToPoint(p));
    poligon_iter++;
  }

  return dist;
}

void Geometry::CreateCube(double x, double width, double length, double H) {
  Point p[8];
  for (int i = 0; i < 4; i++) {
    p[i].x = x;
    p[i + 4].x = x + length;
  }

  p[0].y = p[3].y = p[4].y = p[7].y = -H / 2;
  p[1].y = p[2].y = p[5].y = p[6].y = H / 2;

  p[0].z = p[1].z = p[4].z = p[5].z = width / 2;
  p[2].z = p[3].z = p[6].z = p[7].z = -width / 2;

  auto add_polygon = [&](const Point& p1, const Point& p2, const Point& p3,
                         const Point& normal) {
    auto polygon = std::make_unique<Polygon>(p1, p2, p3, normal);
    polygon->flux = 0;
    polygon->force = Point(0, 0, 0);
    poligons.push_back(std::move(polygon));
  };

  add_polygon(p[0], p[1], p[2], Point(-1, 0, 0));
  add_polygon(p[0], p[2], p[3], Point(-1, 0, 0));
  add_polygon(p[4], p[5], p[6], Point(1, 0, 0));
  add_polygon(p[4], p[6], p[7], Point(1, 0, 0));

  add_polygon(p[1], p[2], p[6], Point(0, 1, 0));
  add_polygon(p[1], p[6], p[5], Point(0, 1, 0));
  add_polygon(p[0], p[3], p[7], Point(0, -1, 0));
  add_polygon(p[0], p[7], p[4], Point(0, -1, 0));

  add_polygon(p[0], p[1], p[5], Point(0, 0, -1));
  add_polygon(p[0], p[5], p[4], Point(0, 0, -1));
  add_polygon(p[3], p[2], p[6], Point(0, 0, 1));
  add_polygon(p[3], p[6], p[7], Point(0, 0, 1));
}

void Geometry::CreateCylinder(double x, double radius, double length,
                              int segments) {
  const double Pi = 3.14159265358979323846;

  // Generate points for the front and back circles
  std::vector<Point> front_circle(segments);
  std::vector<Point> back_circle(segments);

  for (int i = 0; i < segments; i++) {
    double angle = 2.0 * Pi * i / segments;
    double y = radius * cos(angle);
    double z = radius * sin(angle);

    front_circle[i].Set(x, y, z);
    back_circle[i].Set(x + length, y, z);
  }

  // Center points for the circular faces
  Point front_center(x, 0, 0);
  Point back_center(x + length, 0, 0);

  auto add_polygon = [&](const Point& p1, const Point& p2, const Point& p3,
                         const Point& normal) {
    auto polygon = std::make_unique<Polygon>(p1, p2, p3, normal);
    polygon->flux = 0;
    polygon->force = Point(0, 0, 0);
    poligons.push_back(std::move(polygon));
  };

  // Create front circular face (triangles from center to edge)
  for (int i = 0; i < segments; i++) {
    int next_i = (i + 1) % segments;
    Point normal = Point(-1, 0, 0);  // Normal pointing inward
    add_polygon(front_center, front_circle[next_i], front_circle[i], normal);
  }

  // Create back circular face (triangles from center to edge)
  for (int i = 0; i < segments; i++) {
    int next_i = (i + 1) % segments;
    Point normal = Point(1, 0, 0);  // Normal pointing inward
    add_polygon(back_center, back_circle[i], back_circle[next_i], normal);
  }

  // Create cylindrical surface (rectangular segments)
  for (int i = 0; i < segments; i++) {
    int next_i = (i + 1) % segments;

    // Each rectangle is divided into two triangles
    Point p1 = front_circle[i];
    Point p2 = front_circle[next_i];
    Point p3 = back_circle[next_i];
    Point p4 = back_circle[i];

    // Calculate normal for the cylindrical surface (pointing inward)
    Point center_to_edge = p1 - front_center;
    center_to_edge.Normalize();
    Point normal = center_to_edge;

    // First triangle
    add_polygon(p1, p2, p3, normal);
    // Second triangle
    add_polygon(p1, p3, p4, normal);
  }
}

int Geometry::FixPolygons() {
  int res = 0;
  for (auto& polygon : poligons) {
    if (polygon->Fix()) res++;
  }
  return res;
}

void Geometry::ReverseNormals() {
  for (auto& polygon : poligons) {
    polygon->normal = polygon->normal * -1.;
  }
}
}  // namespace mc3d
