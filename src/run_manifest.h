#pragma once
#include <filesystem>
#include <map>
#include <set>

#include "../third_party/json.hpp"
#include "config.h"
#include "geometry.h"
namespace mc3d {
nlohmann::json EffectiveConfig(const SimulationConfig& config);
class RunManifest {
 public:
  RunManifest(SimulationConfig& config, int argc, char** argv);
  void GeometryHash(const Geometry& geometry);
  void RuntimePolicy(double time, std::size_t step, const std::string& mode,
                     std::size_t workers);
  void Progress(double time, std::size_t steps);
  void Finish(const std::string& status, const std::string& error = "");

 private:
  nlohmann::json data_;
  std::filesystem::path directory_, original_cwd_, output_;
  bool isolated_ = false;
  std::map<std::filesystem::path, std::filesystem::file_time_type> previous_;
  std::set<std::filesystem::path> explicit_outputs_;
  void Save();
};
}  // namespace mc3d
