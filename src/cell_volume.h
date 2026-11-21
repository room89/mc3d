#pragma once
#include <array>
#include <cstddef>
#include <optional>

#include "point.h"
namespace mc3d {
struct CellVolumeSummary {
  std::array<std::size_t, 3> counts{};
  double full = 0, gas = 0, solid = 0;
};
enum class CellRegion { Gas, Cut, Solid };
struct CellVolume {
  double full = 0, gas = 0, solid = 0;
  CellRegion region = CellRegion::Gas;
  std::optional<Point> gas_center;
};
}  // namespace mc3d
