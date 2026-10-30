#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "particle.h"

namespace {
using mc3d::Particle;
using mc3d::Point;
template <class T>
constexpr bool ValueType =
    std::is_standard_layout_v<T> && std::is_trivially_copyable_v<T> &&
    std::is_trivially_destructible_v<T> &&
    std::is_trivially_copy_constructible_v<T> &&
    std::is_trivially_copy_assignable_v<T> &&
    std::is_trivially_move_constructible_v<T> &&
    std::is_trivially_move_assignable_v<T>;
static_assert(ValueType<Point> && ValueType<Particle>);
// The deliberately user-provided Point constructor must not become '= default'.
static_assert(!std::is_trivially_default_constructible_v<Point>);
static_assert(sizeof(Point) == 3 * sizeof(double));
static_assert(offsetof(Point, x) == 0 && offsetof(Point, y) == sizeof(double) &&
              offsetof(Point, z) == 2 * sizeof(double));
static_assert(sizeof(Particle) == 2 * sizeof(Point));
static_assert(offsetof(Particle, velocity) == 0 &&
              offsetof(Particle, position) == sizeof(Point));
static_assert(
    std::is_same_v<decltype(std::declval<Point&>() += Point(0, 0, 0)), Point>);
static_assert(
    std::is_same_v<decltype(std::declval<Point&>() -= Point(0, 0, 0)), Point>);
static_assert(std::is_same_v<decltype(std::declval<Point&>() *= 1.), Point>);
static_assert(std::is_same_v<decltype(std::declval<Point&>() /= 1.), Point>);

void Exact(Point a, Point b) {
  EXPECT_EQ(a.x, b.x);
  EXPECT_EQ(a.y, b.y);
  EXPECT_EQ(a.z, b.z);
}
void NearReference(double actual, long double expected, long double scale) {
  const long double tolerance = 64 * std::numeric_limits<double>::epsilon() *
                                    std::max(std::abs(expected), scale) +
                                8 * std::numeric_limits<double>::denorm_min();
  EXPECT_LE(std::abs(static_cast<long double>(actual) - expected), tolerance);
}
}  // namespace

TEST(ValueTypesTest, ConstructorsAccessorsAndVectorRelocationPreserveValues) {
  Particle zero;
  Exact(zero.position, {0, 0, 0});
  Exact(zero.velocity, {0, 0, 0});
  Particle p(1, 2, 3, -4, -5, -6);
  Exact(p.GetPosition(), {1, 2, 3});
  Exact(p.GetVelocity(), {-4, -5, -6});
  Particle q({1, 2, 3}, {-4, -5, -6});
  auto copy = q;
  auto moved = std::move(copy);
  p = moved;
  q = std::move(p);
  Exact(q.position, {1, 2, 3});
  Exact(q.velocity, {-4, -5, -6});
  q.SetX(7);
  q.SetY(8);
  q.SetZ(9);
  q.SetU(-1);
  q.SetV(-2);
  q.SetW(-3);
  EXPECT_EQ(q.GetX(), 7);
  EXPECT_EQ(q.GetY(), 8);
  EXPECT_EQ(q.GetZ(), 9);
  EXPECT_EQ(q.GetU(), -1);
  EXPECT_EQ(q.GetV(), -2);
  EXPECT_EQ(q.GetW(), -3);
  q.SetPosition({1, 2, 3});
  q.SetVelocity({4, 5, 6});
  q.Move(.5);
  Exact(q.position, {3, 4.5, 6});
  std::vector<Particle> values;
  for (int i = 0; i < 1025; ++i)
    values.emplace_back(Point(i, -i, i + .5), Point(-i, i + .25, 2 * i));
  auto relocated = values;
  relocated.reserve(4096);
  auto transferred = std::move(relocated);
  for (size_t i = 0; i < values.size(); ++i) {
    Exact(values[i].position, transferred[i].position);
    Exact(values[i].velocity, transferred[i].velocity);
  }
}

TEST(ValueTypesTest, VectorArithmeticMatchesScaledLongDoubleReference) {
  const std::array<Point, 5> inputs{{{3, -4, 5},
                                     {1e16, 1, -1e16},
                                     {1e100, -2e100, 3e100},
                                     {1e-100, -2e-100, 3e-100},
                                     {-7.5, .125, -.25}}};
  for (auto a : inputs)
    for (auto b : inputs) {
      const std::array<long double, 3> x{a.x, a.y, a.z}, y{b.x, b.y, b.z};
      const auto sum = a + b, difference = a - b, product = a * (-3.25),
                 reversed = (-3.25) * a, quotient = a / (-3.25);
      const std::array<double, 3> sums{sum.x, sum.y, sum.z},
          differences{difference.x, difference.y, difference.z};
      const std::array<double, 3> products{product.x, product.y, product.z},
          quotients{quotient.x, quotient.y, quotient.z};
      Exact(product, reversed);
      long double dot = 0, dot_scale = 0, square = 0;
      for (size_t i = 0; i < 3; ++i) {
        NearReference(sums[i], x[i] + y[i], std::abs(x[i]) + std::abs(y[i]));
        NearReference(differences[i], x[i] - y[i],
                      std::abs(x[i]) + std::abs(y[i]));
        NearReference(products[i], x[i] * (-3.25L), std::abs(x[i] * 3.25L));
        NearReference(quotients[i], x[i] / (-3.25L), std::abs(x[i] / 3.25L));
        dot += x[i] * y[i];
        dot_scale += std::abs(x[i] * y[i]);
        square += x[i] * x[i];
      }
      NearReference(a * b, dot, dot_scale);
      NearReference(a.Mod(), std::sqrt(square), std::sqrt(square));
      NearReference(a.Volume(), x[0] * x[1] * x[2],
                    std::abs(x[0] * x[1] * x[2]));
      auto cross = a.Cross(b);
      const std::array<double, 3> components{cross.x, cross.y, cross.z};
      for (size_t i = 0; i < 3; ++i) {
        const auto j = (i + 1) % 3, k = (i + 2) % 3;
        NearReference(components[i], x[j] * y[k] - x[k] * y[j],
                      std::abs(x[j] * y[k]) + std::abs(x[k] * y[j]));
      }
      auto normalized = a;
      EXPECT_EQ(&normalized.Normalize(), &normalized);
      NearReference(normalized.x, x[0] / std::sqrt(square), 1);
      NearReference(normalized.y, x[1] / std::sqrt(square), 1);
      NearReference(normalized.z, x[2] / std::sqrt(square), 1);
      auto compound = a;
      auto returned = (compound += b);
      Exact(compound, sum);
      Exact(returned, sum);
      returned.x = 0;
      Exact(compound, sum);
      compound = a;
      Exact(compound -= b, difference);
      Exact(compound, difference);
      compound = a;
      Exact(compound *= (-3.25), product);
      Exact(compound, product);
      compound = a;
      Exact(compound /= (-3.25), quotient);
      Exact(compound, quotient);
      Exact(a.Get(), a);
      EXPECT_EQ(a.GetX(), a.x);
      Exact(compound.Set(a.x, a.y, a.z), a);
      Exact(compound, a);
    }
}

TEST(ValueTypesTest, ExceptionalFloatingValuesKeepExistingSemantics) {
  Point p(-0., 0., -0.);
  EXPECT_TRUE(std::signbit(p.x));
  EXPECT_FALSE(std::signbit(p.y));
  EXPECT_TRUE(std::signbit(p.z));
  auto divided = Point(1, -1, 0) / 0.;
  EXPECT_TRUE(std::isinf(divided.x));
  EXPECT_FALSE(std::signbit(divided.x));
  EXPECT_TRUE(std::isinf(divided.y));
  EXPECT_TRUE(std::signbit(divided.y));
  EXPECT_TRUE(std::isnan(divided.z));
  Point zero(0, 0, 0);
  zero.Normalize();
  EXPECT_TRUE(std::isnan(zero.x));
  EXPECT_TRUE(std::isnan(zero.y));
  EXPECT_TRUE(std::isnan(zero.z));
  Exact(!Point(-1, 1e-7, std::nextafter(1e-7, 0.)), {1, 0, 1});
}
