#pragma once

#include <array>
#include <cstddef>
#include <optional>
#include <string>

namespace mc3d {

enum class BoundaryType {
  None = 0,
  Mirror,
  Periodic,
  Free,
  HyperFree,
};

enum class BoundaryFace : std::size_t {
  XNeg = 0,
  XPos,
  YNeg,
  YPos,
  ZNeg,
  ZPos,
};

struct SimulationConfig {
  double Lx = 1.0;
  double Ly = 1.0;
  double Lz = 1.0;

  double Kn = 0.05;
  double Cu = 0.7;
  double temperature = 1.0;
  double alpha = 0.0;

  unsigned int cells_x = 20;
  unsigned int cells_y = 20;
  unsigned int cells_z = 3;
  unsigned int particles_per_cell = 50;

  double S = 10.0;
  double end_time = 1.0;

  std::optional<double> apex_x;
  std::optional<double> apex_y;
  std::optional<double> apex_z;

  bool write_default_snapshot_before_compute = true;
  bool write_default_snapshot_after_compute = true;
  bool write_speed_before = true;
  bool write_speed_after = true;
  bool write_times = true;
  bool snapshots_binary = false;
  double snapshot_interval = 0.0;

  std::string precompute_snapshot = "data_NU.dat";
  std::string final_cell_snapshot = "end_cell.net";
  std::string final_snapshot;

  std::string geometry_type = "none";
  std::string geometry_file;
  bool geometry_fix_polygons = false;
  bool geometry_reverse_normals = false;
  std::optional<double> geometry_fragment_length;
  std::optional<double> geometry_scale;
  std::optional<double> geometry_move_x;
  std::optional<double> geometry_move_y;
  std::optional<double> geometry_move_z;

  double geometry_wedge_x = -0.25;
  double geometry_wedge_width = 0.5;
  double geometry_wedge_length = 0.5;
  double geometry_wedge_alpha = 3.14159265358979323846 / 18.0;

  double geometry_pyramid_x = -0.25;
  double geometry_pyramid_width = 0.5;
  double geometry_pyramid_length = 0.5;
  double geometry_pyramid_height = 0.16;

  double geometry_cube_x = -0.25;
  double geometry_cube_width = 0.5;
  double geometry_cube_length = 0.5;
  double geometry_cube_height = 0.5;

  std::array<BoundaryType, 6> boundary_types{
      BoundaryType::HyperFree, BoundaryType::HyperFree,
      BoundaryType::HyperFree, BoundaryType::HyperFree,
      BoundaryType::HyperFree, BoundaryType::HyperFree};

  std::string config_path;
};

SimulationConfig LoadSimulationConfig(int argc, char* argv[]);

const char* BoundaryTypeToString(BoundaryType type);
const char* BoundaryFaceToString(BoundaryFace face);

constexpr std::size_t ToIndex(BoundaryFace face) {
  return static_cast<std::size_t>(face);
}

}  // namespace mc3d
