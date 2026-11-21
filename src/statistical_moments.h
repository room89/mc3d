#pragma once
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "point.h"
namespace mc3d {
struct ConservedQuantity {
  std::int64_t count = 0;
  Point momentum{0, 0, 0};
  double energy = 0, absolute_momentum = 0;
  void Add(Point v) {
    const double e = .5 * (v * v);
    if (!std::isfinite(e))
      throw std::runtime_error("Nonfinite particle energy in balance");
    ++count;
    momentum += v;
    energy += e;
    absolute_momentum += v.Mod();
  }
  ConservedQuantity& operator+=(const ConservedQuantity& b) {
    count += b.count;
    momentum += b.momentum;
    energy += b.energy;
    absolute_momentum += b.absolute_momentum;
    return *this;
  }
};
struct WindowMoments {
  long double duration = 0, a0 = 0;
  std::array<long double, 3> origin{}, a1{}, a3{};
  std::array<long double, 9> a2{};
  bool initialized = false;
  void Add(Point v, double weight) {
    if (!initialized) {
      origin = {v.x, v.y, v.z};
      initialized = true;
    }
    std::array<long double, 3> c{v.x - origin[0], v.y - origin[1],
                                 v.z - origin[2]};
    long double square = 0;
    for (auto x : c) square += x * x;
    a0 += weight;
    for (int i = 0; i < 3; ++i) {
      a1[i] += weight * c[i];
      a3[i] += weight * c[i] * square;
      for (int j = 0; j < 3; ++j) a2[3 * i + j] += weight * c[i] * c[j];
    }
  }
  struct Result {
    Point velocity{0, 0, 0}, heat{0, 0, 0};
    double temperature = 0;
    std::array<double, 9> pressure{};
  };
  Result Evaluate(double particle_weight, double volume) const {
    Result r;
    if (a0 == 0 || duration == 0 || volume == 0) return r;
    std::array<long double, 3> u{};
    long double tr = 0, u2 = 0;
    for (int i = 0; i < 3; ++i) {
      u[i] = a1[i] / a0;
      tr += a2[3 * i + i];
      u2 += u[i] * u[i];
    }
    r.velocity = {double(u[0] + origin[0]), double(u[1] + origin[1]),
                  double(u[2] + origin[2])};
    r.temperature = double((tr / a0 - u2) / 3);
    const long double scale = particle_weight / (volume * duration);
    std::array<double, 3> q{};
    for (int i = 0; i < 3; ++i) {
      long double second = 0;
      for (int j = 0; j < 3; ++j) {
        second += a2[3 * i + j] * u[j];
        r.pressure[3 * i + j] =
            double(scale * (a2[3 * i + j] - a0 * u[i] * u[j]));
      }
      q[i] = double(.5L * scale *
                    (a3[i] - u[i] * tr - 2 * second + 2 * a0 * u2 * u[i]));
    }
    r.heat = {q[0], q[1], q[2]};
    return r;
  }
};
}  // namespace mc3d
