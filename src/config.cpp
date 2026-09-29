#include "config.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace mc3d {
namespace {

struct OptionContext {
  std::string source;
  std::size_t line = 0;
};

std::string TrimCopy(std::string_view value) {
  const auto begin = value.find_first_not_of(" \t\r\n");
  if (begin == std::string_view::npos) {
    return {};
  }
  const auto end = value.find_last_not_of(" \t\r\n");
  return std::string{value.substr(begin, end - begin + 1)};
}

std::string StripQuotes(std::string_view value) {
  if (value.size() >= 2 && ((value.front() == '"' && value.back() == '"') ||
                            (value.front() == '\'' && value.back() == '\''))) {
    return std::string{value.substr(1, value.size() - 2)};
  }
  return std::string{value};
}

std::string NormalizeKey(std::string_view key) {
  std::string normalized;
  normalized.reserve(key.size());
  for (char ch : key) {
    if (ch == '-' || ch == '.' || ch == ' ' || ch == '/') {
      normalized.push_back('_');
    } else {
      normalized.push_back(
          static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
    }
  }
  return normalized;
}

std::string ToLower(std::string_view value) {
  std::string result;
  result.reserve(value.size());
  for (char ch : value) {
    result.push_back(
        static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  }
  return result;
}

bool IsNoneLike(std::string_view value) {
  const auto lowered = ToLower(TrimCopy(value));
  return lowered == "none" || lowered == "null" || lowered == "off" ||
         lowered == "disable" || lowered == "disabled";
}

bool LooksLikeOptionToken(const std::string& token) {
  if (token.size() >= 2 && token.rfind("--", 0) == 0) {
    return true;
  }
  if (token.empty() || token[0] != '-') {
    return false;
  }
  if (token.size() == 1) {
    return false;
  }
  const unsigned char second = static_cast<unsigned char>(token[1]);
  if (std::isdigit(second) || token[1] == '.') {
    return false;
  }
  return true;
}

std::string_view DescribeOrigin(const OptionContext& ctx, std::string& buffer) {
  if (ctx.source.empty()) {
    return {};
  }
  std::ostringstream oss;
  oss << ctx.source;
  if (ctx.line != 0) {
    oss << ':' << ctx.line;
  }
  buffer = oss.str();
  return buffer;
}

[[noreturn]] void ThrowOptionError(const std::string& key,
                                   const OptionContext& ctx,
                                   const std::string& message) {
  std::string origin_buffer;
  const auto origin = DescribeOrigin(ctx, origin_buffer);

  std::ostringstream oss;
  oss << "Option '" << key << "'";
  if (!origin.empty()) {
    oss << " at " << origin;
  }
  oss << ": " << message;
  throw std::runtime_error(oss.str());
}

double ParseDouble(const std::string& value, const std::string& key,
                   const OptionContext& ctx) {
  try {
    size_t consumed = 0;
    const double result = std::stod(value, &consumed);
    if (consumed != value.size()) {
      ThrowOptionError(key, ctx,
                       "unexpected trailing characters in '" + value + "'");
    }
    return result;
  } catch (const std::exception&) {
    ThrowOptionError(key, ctx,
                     "expected floating-point number, got '" + value + "'");
  }
}

unsigned int ParseUnsigned(const std::string& value, const std::string& key,
                           const OptionContext& ctx) {
  try {
    size_t consumed = 0;
    const auto parsed = std::stoul(value, &consumed);
    if (consumed != value.size()) {
      ThrowOptionError(key, ctx,
                       "unexpected trailing characters in '" + value + "'");
    }
    return static_cast<unsigned int>(parsed);
  } catch (const std::exception&) {
    ThrowOptionError(key, ctx,
                     "expected positive integer value, got '" + value + "'");
  }
}

bool ParseBool(const std::string& value, const std::string& key,
               const OptionContext& ctx) {
  const auto lowered = ToLower(value);
  if (lowered == "true" || lowered == "1" || lowered == "yes" ||
      lowered == "on") {
    return true;
  }
  if (lowered == "false" || lowered == "0" || lowered == "no" ||
      lowered == "off") {
    return false;
  }
  ThrowOptionError(key, ctx, "expected boolean value, got '" + value + "'");
}

BoundaryType ParseBoundaryType(const std::string& value, const std::string& key,
                               const OptionContext& ctx) {
  const auto normalized = NormalizeKey(value);
  if (normalized == "none") {
    return BoundaryType::None;
  }
  if (normalized == "mirror") {
    return BoundaryType::Mirror;
  }
  if (normalized == "periodic") {
    return BoundaryType::Periodic;
  }
  if (normalized == "free") {
    return BoundaryType::Free;
  }
  if (normalized == "hyperfree" || normalized == "hyper_free" ||
      normalized == "hyper-free" || normalized == "giperfree" ||
      normalized == "giper_free") {
    return BoundaryType::HyperFree;
  }
  ThrowOptionError(
      key, ctx,
      "unsupported boundary type '" + value +
          "'. Use one of: none, mirror, periodic, free, hyperfree.");
}

std::string ParseStringValue(const std::string& value) {
  auto trimmed = TrimCopy(value);
  trimmed = StripQuotes(trimmed);
  if (IsNoneLike(trimmed)) {
    return {};
  }
  return trimmed;
}

void ApplyOption(SimulationConfig& config, const std::string& raw_key,
                 const std::string& raw_value, const OptionContext& ctx) {
  const auto key = NormalizeKey(raw_key);
  const auto prepared_value = StripQuotes(TrimCopy(raw_value));

  if (key == "lx") {
    config.Lx = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "ly") {
    config.Ly = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "lz") {
    config.Lz = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "kn") {
    config.Kn = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "cu") {
    config.Cu = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "temperature" || key == "t") {
    config.temperature = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "wall_temperature") {
    config.wall_temperature = ParseDouble(prepared_value, raw_key, ctx);
    if (!(config.wall_temperature > 0) || !std::isfinite(config.wall_temperature))
      throw std::invalid_argument("Wall temperature must be finite and positive");
  } else if (key == "alpha") {
    config.alpha = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "ncx" || key == "cells_x") {
    config.cells_x = ParseUnsigned(prepared_value, raw_key, ctx);
  } else if (key == "ncy" || key == "cells_y") {
    config.cells_y = ParseUnsigned(prepared_value, raw_key, ctx);
  } else if (key == "ncz" || key == "cells_z") {
    config.cells_z = ParseUnsigned(prepared_value, raw_key, ctx);
  } else if (key == "np" || key == "particles_per_cell" || key == "particles") {
    config.particles_per_cell = ParseUnsigned(prepared_value, raw_key, ctx);
  } else if (key == "thread_pool_size" || key == "thread_count" ||
             key == "threads") {
    const auto parsed = ParseUnsigned(prepared_value, raw_key, ctx);
    if (parsed == 0) {
      ThrowOptionError(raw_key, ctx, "thread pool size must be positive");
    }
    config.thread_pool_size = parsed;
  } else if (key == "s" || key == "sigma") {
    config.S = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "end_time" || key == "t_end") {
    config.end_time = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "apex_x") {
    config.apex_x = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "apex_y") {
    config.apex_y = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "apex_z") {
    config.apex_z = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "precompute_snapshot" || key == "initial_snapshot") {
    config.precompute_snapshot = ParseStringValue(prepared_value);
  } else if (key == "final_cell_snapshot") {
    config.final_cell_snapshot = ParseStringValue(prepared_value);
  } else if (key == "final_snapshot") {
    config.final_snapshot = ParseStringValue(prepared_value);
  } else if (key == "write_default_snapshot_before_compute") {
    config.write_default_snapshot_before_compute =
        ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "write_default_snapshot_after_compute") {
    config.write_default_snapshot_after_compute =
        ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "write_speed_before") {
    config.write_speed_before = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "write_speed_after") {
    config.write_speed_after = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "write_times") {
    config.write_times = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "snapshots_binary" || key == "snapshot_binary" ||
             key == "write_binary_snapshot") {
    config.snapshots_binary = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "snapshot_interval" || key == "snapshots_interval" ||
             key == "snapshot_time_interval") {
    config.snapshot_interval = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_type") {
    config.geometry_type = NormalizeKey(prepared_value);
  } else if (key == "geometry_file") {
    config.geometry_file = ParseStringValue(prepared_value);
  } else if (key == "geometry_fix_polygons") {
    config.geometry_fix_polygons = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "geometry_reverse_normals") {
    config.geometry_reverse_normals = ParseBool(prepared_value, raw_key, ctx);
  } else if (key == "geometry_fragment_length") {
    config.geometry_fragment_length = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_scale") {
    config.geometry_scale = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_move_x") {
    config.geometry_move_x = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_move_y") {
    config.geometry_move_y = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_move_z") {
    config.geometry_move_z = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_wedge_x") {
    config.geometry_wedge_x = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_wedge_width") {
    config.geometry_wedge_width = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_wedge_length") {
    config.geometry_wedge_length = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_wedge_alpha") {
    config.geometry_wedge_alpha = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_pyramid_x") {
    config.geometry_pyramid_x = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_pyramid_width") {
    config.geometry_pyramid_width = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_pyramid_length") {
    config.geometry_pyramid_length = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_pyramid_height") {
    config.geometry_pyramid_height = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_cube_x") {
    config.geometry_cube_x = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_cube_width") {
    config.geometry_cube_width = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_cube_length") {
    config.geometry_cube_length = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "geometry_cube_height") {
    config.geometry_cube_height = ParseDouble(prepared_value, raw_key, ctx);
  } else if (key == "boundary_x_neg") {
    config.boundary_types[ToIndex(BoundaryFace::XNeg)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else if (key == "boundary_x_pos" || key == "boundary_x_plus") {
    config.boundary_types[ToIndex(BoundaryFace::XPos)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else if (key == "boundary_y_neg") {
    config.boundary_types[ToIndex(BoundaryFace::YNeg)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else if (key == "boundary_y_pos" || key == "boundary_y_plus") {
    config.boundary_types[ToIndex(BoundaryFace::YPos)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else if (key == "boundary_z_neg") {
    config.boundary_types[ToIndex(BoundaryFace::ZNeg)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else if (key == "boundary_z_pos" || key == "boundary_z_plus") {
    config.boundary_types[ToIndex(BoundaryFace::ZPos)] =
        ParseBoundaryType(prepared_value, raw_key, ctx);
  } else {
    ThrowOptionError(raw_key, ctx, "unknown option");
  }
}

std::string RemoveTrailingComments(std::string_view line) {
  bool in_single_quote = false;
  bool in_double_quote = false;

  for (std::size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];
    if (ch == '\'' && !in_double_quote) {
      in_single_quote = !in_single_quote;
    } else if (ch == '"' && !in_single_quote) {
      in_double_quote = !in_double_quote;
    }

    if (in_single_quote || in_double_quote) {
      continue;
    }

    if (ch == '#') {
      return std::string{line.substr(0, i)};
    }

    if (ch == '/' && i + 1 < line.size() && line[i + 1] == '/') {
      return std::string{line.substr(0, i)};
    }
  }

  return std::string{line};
}

void LoadKeyValueConfig(SimulationConfig& config,
                        const std::filesystem::path& path, bool allow_colon) {
  std::ifstream stream(path);
  if (!stream) {
    throw std::runtime_error("Failed to open config file '" + path.string() +
                             "'");
  }

  std::string line;
  std::size_t line_number = 0;

  while (std::getline(stream, line)) {
    ++line_number;

    std::string prepared = RemoveTrailingComments(line);
    prepared = TrimCopy(prepared);
    if (prepared.empty()) {
      continue;
    }

    std::string key;
    std::string value;

    const auto eq_pos = prepared.find('=');
    std::size_t colon_pos = std::string::npos;
    if (allow_colon) {
      colon_pos = prepared.find(':');
      if (colon_pos != std::string::npos && eq_pos != std::string::npos &&
          colon_pos > eq_pos) {
        colon_pos = std::string::npos;
      }
    }

    std::size_t delimiter_pos = std::string::npos;
    if (eq_pos != std::string::npos &&
        (colon_pos == std::string::npos || eq_pos < colon_pos)) {
      delimiter_pos = eq_pos;
    } else if (colon_pos != std::string::npos) {
      delimiter_pos = colon_pos;
    }

    if (delimiter_pos != std::string::npos) {
      key = TrimCopy(prepared.substr(0, delimiter_pos));
      value = TrimCopy(prepared.substr(delimiter_pos + 1));
    } else {
      std::istringstream iss(prepared);
      if (!(iss >> key)) {
        continue;
      }
      if (!(iss >> value)) {
        OptionContext ctx{path.string(), line_number};
        ThrowOptionError(
            key, ctx,
            "missing value (use 'key value', 'key = value' or 'key: value')");
      }
      std::string rest;
      if (std::getline(iss, rest)) {
        const auto trimmed_rest = TrimCopy(rest);
        if (!trimmed_rest.empty()) {
          value += ' ' + trimmed_rest;
        }
      }
    }

    OptionContext ctx{path.string(), line_number};
    ApplyOption(config, key, value, ctx);
  }
}

class JsonReader {
 public:
  explicit JsonReader(std::string data, std::string source)
      : data_(std::move(data)), source_(std::move(source)) {}

  void ParseObject(const std::function<void(std::string key, std::string value,
                                            const OptionContext&)>& on_entry) {
    SkipWhitespace();
    Expect('{');
    SkipWhitespace();
    if (Match('}')) {
      return;
    }

    while (true) {
      SkipWhitespace();
      const auto key = ParseString();
      SkipWhitespace();
      Expect(':');
      SkipWhitespace();
      const auto value = ParseValue();

      OptionContext ctx{source_, CurrentLine()};
      on_entry(key, value, ctx);

      SkipWhitespace();
      if (Match('}')) {
        break;
      }
      Expect(',');
    }
  }

 private:
  bool Match(char expected) {
    if (pos_ < data_.size() && data_[pos_] == expected) {
      ++pos_;
      return true;
    }
    return false;
  }

  void Expect(char expected) {
    if (!Match(expected)) {
      Throw("expected '" + std::string(1, expected) + "'");
    }
  }

  void SkipWhitespace() {
    while (pos_ < data_.size()) {
      const char ch = data_[pos_];
      if (ch == ' ' || ch == '\t' || ch == '\n' || ch == '\r') {
        ++pos_;
      } else {
        break;
      }
    }
  }

  static bool IsDigit(char ch) { return ch >= '0' && ch <= '9'; }

  std::string ParseString() {
    Expect('"');
    std::string result;
    while (pos_ < data_.size()) {
      char ch = data_[pos_++];
      if (ch == '"') {
        return result;
      }
      if (ch == '\\') {
        if (pos_ >= data_.size()) {
          Throw("unterminated escape sequence");
        }
        char esc = data_[pos_++];
        switch (esc) {
          case '"':
          case '\\':
          case '/':
            result.push_back(esc);
            break;
          case 'b':
            result.push_back('\b');
            break;
          case 'f':
            result.push_back('\f');
            break;
          case 'n':
            result.push_back('\n');
            break;
          case 'r':
            result.push_back('\r');
            break;
          case 't':
            result.push_back('\t');
            break;
          case 'u':
            Throw("unicode escapes are not supported");
          default:
            Throw("invalid escape character");
        }
      } else {
        result.push_back(ch);
      }
    }
    Throw("unterminated string literal");
  }

  std::string ParseNumber() {
    const std::size_t start = pos_;
    if (Match('-')) {
      if (pos_ >= data_.size() || !IsDigit(data_[pos_])) {
        Throw("invalid number");
      }
    }

    if (Match('0')) {
      if (pos_ < data_.size() && IsDigit(data_[pos_])) {
        Throw("leading zeros are not allowed");
      }
    } else {
      while (pos_ < data_.size() && IsDigit(data_[pos_])) {
        ++pos_;
      }
    }

    if (Match('.')) {
      if (pos_ >= data_.size() || !IsDigit(data_[pos_])) {
        Throw("invalid fractional part");
      }
      while (pos_ < data_.size() && IsDigit(data_[pos_])) {
        ++pos_;
      }
    }

    if (Match('e') || Match('E')) {
      if (!Match('+') && !Match('-')) {
        --pos_;
      }
      if (pos_ >= data_.size() || !IsDigit(data_[pos_])) {
        Throw("invalid exponent");
      }
      while (pos_ < data_.size() && IsDigit(data_[pos_])) {
        ++pos_;
      }
    }

    return data_.substr(start, pos_ - start);
  }

  std::string ParseLiteral(std::string_view literal, std::string value) {
    if (data_.compare(pos_, literal.size(), literal) != 0) {
      Throw("expected literal '" + std::string(literal) + "'");
    }
    pos_ += literal.size();
    return value;
  }

  std::string ParseValue() {
    if (pos_ >= data_.size()) {
      Throw("unexpected end of input");
    }
    const char ch = data_[pos_];
    if (ch == '"') {
      return ParseString();
    }
    if (ch == '-' || IsDigit(ch)) {
      return ParseNumber();
    }
    if (ch == 't') {
      return ParseLiteral("true", "true");
    }
    if (ch == 'f') {
      return ParseLiteral("false", "false");
    }
    if (ch == 'n') {
      return ParseLiteral("null", "null");
    }
    Throw("unsupported JSON value");
  }

  [[noreturn]] void Throw(const std::string& message) const {
    std::ostringstream oss;
    oss << "JSON parse error in '" << source_ << "' at line " << CurrentLine()
        << ", column " << CurrentColumn() << ": " << message;
    throw std::runtime_error(oss.str());
  }

  std::size_t CurrentLine() const {
    std::size_t line = 1;
    for (std::size_t i = 0; i < pos_ && i < data_.size(); ++i) {
      if (data_[i] == '\n') {
        ++line;
      }
    }
    return line;
  }

  std::size_t CurrentColumn() const {
    std::size_t column = 1;
    for (std::size_t i = pos_; i > 0; --i) {
      if (data_[i - 1] == '\n') {
        break;
      }
      ++column;
    }
    return column;
  }

  std::string data_;
  std::string source_;
  std::size_t pos_ = 0;
};

void LoadJsonConfig(SimulationConfig& config,
                    const std::filesystem::path& path) {
  std::ifstream stream(path);
  if (!stream) {
    throw std::runtime_error("Failed to open config file '" + path.string() +
                             "'");
  }
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  JsonReader reader(buffer.str(), path.string());

  reader.ParseObject(
      [&](std::string key, std::string value, const OptionContext& ctx) {
        ApplyOption(config, key, value, ctx);
      });
}

void LoadConfigFile(SimulationConfig& config,
                    const std::filesystem::path& path) {
  const auto extension = path.extension().string();
  const auto extension_lower = ToLower(extension);
  if (extension_lower == ".json") {
    LoadJsonConfig(config, path);
    return;
  }
  if (extension_lower == ".yaml" || extension_lower == ".yml") {
    LoadKeyValueConfig(config, path, /*allow_colon=*/true);
    return;
  }

  LoadKeyValueConfig(config, path, /*allow_colon=*/true);
}

void PrintUsage(const char* program) {
  std::cout
      << "Usage: " << (program ? program : "MC3dSolver")
      << " [--config FILE] [--option value]...\n\n"
      << "Key options:\n"
      << "  --lx, --ly, --lz                 Domain size along axes\n"
      << "  --ncx, --ncy, --ncz              Number of cells along axes\n"
      << "  --np                             Particles per cell\n"
      << "  --kn, --cu, --temperature        Physical parameters\n"
      << "  --wall-temperature VALUE        Diffuse body temperature (default 1)\n"
      << "  --s                              Accommodation coefficient S\n"
      << "  --end-time                       Simulation end time\n"
      << "  --snapshot-interval TIME         Snapshot interval in time units\n"
      << "  --thread-pool-size COUNT         Number of worker threads\n"
      << "  --snapshots-binary               Enable binary snapshot output\n"
      << "  --boundary-x-neg TYPE            Boundary type (none, mirror, "
         "periodic, free, hyperfree)\n"
      << "  --geometry-type TYPE             Body geometry (none, wedge, "
         "pyramid, cube)\n"
      << "  --geometry-file PATH             Load geometry from STL file\n"
      << "\n"
      << "All options can be provided in a configuration file using "
         "'key = value'.\n";
}

}  // namespace

SimulationConfig LoadSimulationConfig(int argc, char* argv[]) {
  SimulationConfig config;
  if (argc <= 1 || argv == nullptr) {
    return config;
  }

  std::vector<std::string> args;
  args.reserve(static_cast<std::size_t>(argc - 1));
  for (int i = 1; i < argc; ++i) {
    args.emplace_back(argv[i]);
  }

  std::vector<std::filesystem::path> config_paths;

  for (std::size_t i = 0; i < args.size(); ++i) {
    const auto& arg = args[i];
    if (arg == "--help" || arg == "-h") {
      PrintUsage(argv[0]);
      std::exit(EXIT_SUCCESS);
    }
    if (arg == "--config" || arg == "-c") {
      if (i + 1 >= args.size()) {
        throw std::runtime_error("Option '--config' requires a path argument");
      }
      config_paths.emplace_back(args[i + 1]);
      ++i;
      continue;
    }
    constexpr std::string_view config_eq = "--config=";
    if (arg.rfind(config_eq, 0) == 0) {
      config_paths.emplace_back(arg.substr(config_eq.size()));
    }
  }

  for (const auto& path : config_paths) {
    LoadConfigFile(config, path);
    config.config_path = path.string();
  }

  for (std::size_t i = 0; i < args.size(); ++i) {
    const auto& arg = args[i];
    if (arg == "--help" || arg == "-h") {
      continue;
    }
    if (arg == "--config" || arg == "-c") {
      ++i;
      continue;
    }
    constexpr std::string_view config_eq = "--config=";
    if (arg.rfind(config_eq, 0) == 0) {
      continue;
    }

    if (arg.size() > 2 && arg.rfind("--", 0) == 0) {
      std::string key;
      std::string value;
      const auto eq_pos = arg.find('=');
      if (eq_pos != std::string::npos) {
        key = arg.substr(2, eq_pos - 2);
        value = arg.substr(eq_pos + 1);
      } else {
        key = arg.substr(2);
        if (i + 1 < args.size() && !LooksLikeOptionToken(args[i + 1])) {
          value = args[++i];
        } else {
          value = "true";
        }
      }

      OptionContext ctx{"command line", 0};
      ApplyOption(config, key, value, ctx);
      continue;
    }

    if (arg.size() > 1 && arg[0] == '-') {
      throw std::runtime_error("Unknown short option '" + arg + "'");
    }

    throw std::runtime_error("Unexpected positional argument '" + arg + "'");
  }

  return config;
}

const char* BoundaryTypeToString(BoundaryType type) {
  switch (type) {
    case BoundaryType::None:
      return "none";
    case BoundaryType::Mirror:
      return "mirror";
    case BoundaryType::Periodic:
      return "periodic";
    case BoundaryType::Free:
      return "free";
    case BoundaryType::HyperFree:
      return "hyperfree";
  }
  return "unknown";
}

const char* BoundaryFaceToString(BoundaryFace face) {
  switch (face) {
    case BoundaryFace::XNeg:
      return "x-";
    case BoundaryFace::XPos:
      return "x+";
    case BoundaryFace::YNeg:
      return "y-";
    case BoundaryFace::YPos:
      return "y+";
    case BoundaryFace::ZNeg:
      return "z-";
    case BoundaryFace::ZPos:
      return "z+";
  }
  return "unknown";
}

}  // namespace mc3d
