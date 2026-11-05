#include "run_manifest.h"

#include <chrono>
#include <fstream>

#include "../third_party/picosha2.h"
#include "build_info.h"
#include "utils/utils.hpp"
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif
namespace mc3d {
namespace fs = std::filesystem;
using Json = nlohmann::json;
namespace {
std::string Hash(const fs::path& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("Cannot hash " + path.string());
  picosha2::hash256_one_by_one hasher;
  char buf[65536];
  while (in.read(buf, sizeof(buf)) || in.gcount())
    hasher.process(buf, buf + in.gcount());
  if (!in.eof()) throw std::runtime_error("Cannot read " + path.string());
  hasher.finish();
  return picosha2::get_hash_hex_string(hasher);
}
void CheckUnits(const Json& units, const SimulationConfig& c) {
  if (units.at("schema_version") != 1)
    throw std::runtime_error("Unsupported units schema");
  const auto& s = units.at("scales");
  for (auto it = s.begin(); it != s.end(); ++it) {
    const double value = it.value().get<double>();
    if (!(value > 0) || !std::isfinite(value))
      throw std::runtime_error("Invalid scale: " + it.key());
  }
  auto check = [](double actual, double expected) {
    if (!std::isfinite(expected) || !(expected > 0) ||
        std::abs(actual / expected - 1) > 1e-12)
      throw std::runtime_error("Inconsistent physical unit scales");
  };
  double l = s.at("length_m"), t = s.at("temperature_k"),
         v = s.at("velocity_m_s"), n = s.at("number_density_m3");
  double m = units.at("reference_gas").at("mass_kg"), rho = m * n;
  check(v, std::sqrt(1.380649e-23 * t / m));
  check(s.at("mass_density_kg_m3"), rho);
  check(s.at("time_s"), l / v);
  check(s.at("area_m2"), l * l);
  check(s.at("pressure_pa"), rho * v * v);
  check(s.at("force_n"), rho * v * v * l * l);
  check(s.at("energy_j"), rho * v * v * l * l * l);
  check(s.at("particle_energy_j"), m * v * v);
  check(s.at("heat_capacity_j_k"), rho * v * v * l * l * l / t);
  double count =
      double(c.particles_per_cell) * c.cells_x * c.cells_y * c.cells_z;
  check(s.at("molecules_per_model_particle"),
        n * c.Lx * c.Ly * c.Lz * l * l * l / count);
  if (c.collisions.kind == CollisionKind::VHS) {
    check(c.collisions.mass, m);
    check(c.collisions.temperature_unit, t);
  }
  if (c.collisions.kind != CollisionKind::None) {
    const auto& gas = units.at("reference_gas");
    double w = gas.at("omega"), tr = gas.at("reference_temperature_k"),
           d = gas.at("diameter_m");
    double mean = std::acos(-1.) * d * d / std::tgamma(2.5 - w) *
                  std::sqrt(4 * 1.380649e-23 * tr / m) * std::tgamma(2.5 - w) /
                  std::tgamma(1.5);
    check(c.Kn, std::sqrt(8 * 1.380649e-23 * tr / (std::acos(-1.) * m)) /
                    (n * mean * c.Lx * l));
  }
}
double Now() {
  return std::chrono::duration<double>(
             std::chrono::system_clock::now().time_since_epoch())
      .count();
}
std::string Executable() {
#ifdef __APPLE__
  uint32_t size = 0;
  _NSGetExecutablePath(nullptr, &size);
  std::string buffer(size, '\0');
  if (_NSGetExecutablePath(buffer.data(), &size))
    throw std::runtime_error("Cannot locate executable");
  return fs::canonical(buffer.c_str()).string();
#else
  return fs::canonical("/proc/self/exe").string();
#endif
}
}  // namespace
Json EffectiveConfig(const SimulationConfig& c) {
  Json j = Json::object();
#include "effective_config.inc"
  for (size_t i = 0; i < 6; ++i) {
    std::string key = "boundary_";
    const char* names[] = {"x_neg", "x_pos", "y_neg",
                           "y_pos", "z_neg", "z_pos"};
    key += names[i];
    j[key] = BoundaryTypeToString(c.boundary_types[i]);
  }
  const auto& p = c.collisions;
  j["collision_model"] = p.kind == CollisionKind::HS    ? "hs"
                         : p.kind == CollisionKind::VHS ? "vhs"
                                                        : "none";
  if (p.kind == CollisionKind::VHS) {
    j["molecular_mass_kg"] = p.mass;
    j["vhs_diameter_m"] = p.diameter;
    j["vhs_reference_temperature_k"] = p.reference_temperature;
    j["vhs_omega"] = p.omega;
    j["temperature_unit_k"] = p.temperature_unit;
  }
  return j;
}
RunManifest::RunManifest(SimulationConfig& c, int argc, char** argv) {
  original_cwd_ = fs::current_path();
  isolated_ = !c.output_dir.empty();
  output_ = isolated_ ? fs::absolute(c.output_dir) : original_cwd_;
  if (isolated_ && fs::exists(output_) && !fs::is_empty(output_))
    throw std::runtime_error("Output directory must be empty");
  if (isolated_)
    for (const auto& name :
         {c.precompute_snapshot, c.final_cell_snapshot, c.final_snapshot}) {
      fs::path p(name);
      auto normalized = p.lexically_normal();
      if (normalized == "manifest.json" || normalized == "manifest.json.tmp" ||
          (!normalized.empty() && (*normalized.begin() == "inputs" ||
                                   *normalized.begin() == "run_metadata")))
        throw std::runtime_error("Output path conflicts with run metadata");
      if (p.is_absolute())
        throw std::runtime_error("Absolute output path forbidden");
      for (const auto& part : p)
        if (part == "..")
          throw std::runtime_error("Output path escapes directory");
    }
  if (!c.geometry_file.empty())
    c.geometry_file = fs::absolute(c.geometry_file).string();
  if (!c.units_file.empty()) c.units_file = fs::absolute(c.units_file).string();
  if (!c.provenance_file.empty())
    c.provenance_file = fs::absolute(c.provenance_file).string();
  const auto id = std::to_string(
      std::chrono::system_clock::now().time_since_epoch().count());
  directory_ = isolated_ ? output_ : output_ / "run_metadata" / id;
  fs::create_directories(directory_ / "inputs");
  for (const auto& name :
       {c.precompute_snapshot, c.final_cell_snapshot, c.final_snapshot})
    if (!name.empty()) {
      explicit_outputs_.insert(fs::absolute(output_ / name));
      explicit_outputs_.insert(fs::absolute(output_ / (name + ".surface.csv")));
    }
  for (auto& entry : fs::directory_iterator(output_))
    if (entry.is_regular_file())
      previous_[entry.path()] = entry.last_write_time();
  for (auto& path : explicit_outputs_)
    if (fs::is_regular_file(path)) previous_[path] = fs::last_write_time(path);
  data_ = {{"schema_version", 1},
           {"status", "running"},
           {"started_unix", Now()},
           {"effective_config", EffectiveConfig(c)},
           {"seed", *c.seed},
           {"thread_pool_size", c.thread_pool_size},
           {"rng", utils::RandomStream::kAlgorithm},
           {"units_mode", "normalized_only"},
           {"build",
            {{"revision", MC3D_REV},
             {"dirty", MC3D_DIRTY},
             {"compiler", MC3D_COMPILER},
             {"type", MC3D_BUILD_TYPE},
             {"flags", MC3D_BUILD_FLAGS},
             {"release_flags", MC3D_RELEASE_FLAGS},
             {"debug_flags", MC3D_DEBUG_FLAGS}}},
           {"executable", Executable()},
           {"binary_sha256", Hash(Executable())},
           {"working_directory", original_cwd_.string()},
           {"output_directory", output_.string()},
           {"argv", Json::array()},
           {"inputs", Json::array()}};
  if (std::string(MC3D_REV) == "unknown") data_["build"]["dirty"] = nullptr;
  for (int i = 0; i < argc; ++i) data_["argv"].push_back(argv[i]);
  auto copy = [&](const std::string& source, const std::string& role) {
    auto destination = directory_ / "inputs" /
                       (std::to_string(data_["inputs"].size()) + "_" +
                        fs::path(source).filename().string());
    fs::copy_file(source, destination);
    data_["inputs"].push_back(
        {{"role", role},
         {"original", fs::absolute(source).string()},
         {"path", fs::relative(destination, directory_).string()},
         {"sha256", Hash(destination)}});
    return destination.string();
  };
  Save();
  try {
    for (auto& path : c.input_configs) copy(path, "config");
    if (!c.geometry_file.empty())
      data_["effective_config"]["geometry_file"] =
          copy(c.geometry_file, "geometry");
    if (!c.units_file.empty()) {
      std::ifstream f(c.units_file);
      Json units;
      f >> units;
      CheckUnits(units, c);
      if (units.contains("source_file"))
        copy(units["source_file"].get<std::string>(), "physical_input");
      auto effective = EffectiveConfig(c);
      for (auto it = units.at("normalized_parameters").begin();
           it != units.at("normalized_parameters").end(); ++it) {
        if (!effective.contains(it.key()) || effective[it.key()] != it.value())
          throw std::runtime_error("Units mismatch: " + it.key());
      }
      data_["provenance"] = units.value("provenance", Json());
      data_["units"] = units;
      data_["units_mode"] = "si";
      copy(c.units_file, "units");
    }
    if (!c.provenance_file.empty()) {
      std::ifstream stream(c.provenance_file);
      Json provenance;
      stream >> provenance;
      data_["provenance"] = provenance;
      copy(c.provenance_file, "provenance");
    }
    if (data_.contains("provenance") && data_["provenance"].is_object() &&
        data_["provenance"].contains("files"))
      for (const auto& file : data_["provenance"]["files"])
        copy(file.get<std::string>(), "reference_data");
    if (!c.geometry_file.empty())
      c.geometry_file =
          data_["effective_config"]["geometry_file"].get<std::string>();
    Save();
    if (isolated_) fs::current_path(output_);
  } catch (const std::exception& error) {
    data_["status"] = "failed";
    data_["error"] = error.what();
    Save();
    throw;
  }
}
void RunManifest::Save() {
  auto temp = directory_ / "manifest.json.tmp";
  std::ofstream f(temp);
  f << data_.dump(2) << '\n';
  f.close();
  if (!f) throw std::runtime_error("Cannot write manifest");
  fs::rename(temp, directory_ / "manifest.json");
}
void RunManifest::GeometryHash(const Geometry& geometry) {
  Json points = Json::array();
  for (size_t i = 0; i < geometry.PolygonCount(); ++i) {
    auto& p = geometry.GetPolygon(i);
    for (auto v : {p.p1, p.p2, p.p3}) points.push_back({v.x, v.y, v.z});
  }
  data_["geometry_sha256"] = picosha2::hash256_hex_string(points.dump());
  Save();
}
void RunManifest::Progress(double t, size_t steps) {
  data_["simulation_time"] = t;
  data_["steps"] = steps;
  Save();
}
void RunManifest::Finish(const std::string& status, const std::string& error) {
  data_["status"] = status;
  data_["error"] = error;
  data_["finished_unix"] = Now();
  data_["outputs"] = Json::array();
  std::set<fs::path> candidates = explicit_outputs_;
  for (auto& entry : fs::directory_iterator(output_))
    if (entry.is_regular_file()) candidates.insert(entry.path());
  for (const auto& path : candidates) {
    if (!fs::is_regular_file(path)) continue;
    if (previous_.contains(path) &&
        previous_.at(path) == fs::last_write_time(path))
      continue;
    auto relative = fs::relative(path, output_);
    auto name = relative.generic_string();
    if (name.starts_with("run_metadata/") || name.starts_with("inputs/") ||
        name == "manifest.json" || name == "manifest.json.tmp")
      continue;
    // Only solver outputs, not unrelated files in a legacy working directory.
    auto ext = path.extension();
    if (!explicit_outputs_.contains(path) && ext != ".dat" && ext != ".csv" &&
        ext != ".net")
      continue;
    data_["outputs"].push_back(
        {{"path", name},
         {"sha256", Hash(path)},
         {"kind", name == "time.dat" ? "timing" : "physical"}});
  }
  Save();
  if (isolated_) fs::current_path(original_cwd_);
}
}  // namespace mc3d
