#include "simulation_runner.h"
#include <charconv>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include "utils/logger.hpp"
#include "utils/stopwatch.hpp"
namespace mc3d {
namespace {

constexpr std::size_t kEstimatedRowSize = 90;
constexpr int kOutputPrecision = 8;
constexpr double kSnapshotTolerance = 1e-9;

template <typename Int>
void AppendIntegral(std::string& out, Int value) {
  char buffer[32];
  auto result = std::to_chars(std::begin(buffer), std::end(buffer), value);
  if (result.ec == std::errc()) {
    out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
  } else {
    out.append(std::to_string(value));
  }
}

void AppendDouble(std::string& out, double value) {
  char buffer[64];
  auto result = std::to_chars(std::begin(buffer), std::end(buffer), value,
                              std::chars_format::general, kOutputPrecision);
  if (result.ec == std::errc()) {
    out.append(buffer, static_cast<std::size_t>(result.ptr - buffer));
  } else {
    std::ostringstream fallback;
    fallback.setf(std::ios::scientific);
    fallback.precision(kOutputPrecision > 0 ? kOutputPrecision - 1 : 0);
    fallback << value;
    out.append(fallback.str());
  }
}

}  // namespace

void SnapshotWriter::SetBinaryOutput(bool enabled) { binary_output_ = enabled; }

bool SnapshotWriter::WriteFile(const std::string& file_name) {
  std::ofstream file(file_name, std::ios::out | std::ios::binary);
  if (!file) {
    LOG_ERROR() << "Failed to open output file: " << file_name;
    return false;
  }

  const bool success = binary_output_ ? WriteBinarySnapshot(file) : WriteTextSnapshot(file);
  file.close();
  if (!success || !file) return false;
  if (cluster_.body_ && cluster_.body_->SurfaceHeatingEnabled())
    return cluster_.body_->WriteSurfaceSnapshot(file_name + ".surface.csv", cluster_.t);
  return true;
}

bool SnapshotWriter::WriteFile() {
  if (cluster_.body_ && cluster_.body_->SurfaceHeatingEnabled() &&
      !cluster_.body_->WriteSurfaceSnapshot("data.dat.surface.csv", cluster_.t)) return false;
  std::ofstream file1("data.dat");
  for (auto& cell : cluster_.cells) {
    auto ap = cell.GetCenter();

    file1 << ap << "\t"
          << double(cell.GetParticleCount()) / (cell.GetVolume() * cluster_.density)
          << "\t" << cell.ComputeMoments().temperature << '\n';
  }
  file1.close();
  return file1.good();
}

bool SnapshotWriter::WriteTextSnapshot(std::ofstream& file) {
  const std::size_t cell_count = cluster_.cells.size();
  std::string buffer;
  buffer.reserve(1 + kEstimatedRowSize * cell_count);
  buffer.append("x;y;z;N;ro;T;vx;vy;vz;E\n");

  for (auto& cell : cluster_.cells) {
    const auto center = cell.GetCenter();
    const auto moments = cell.ComputeMoments();
    const auto velocity = moments.velocity;
    const unsigned int particle_count = cell.GetParticleCount();
    const double density_value =
        cell.GetVolume() > 0 && cluster_.density != 0.0
            ? static_cast<double>(particle_count) / (cell.GetVolume() * cluster_.density)
            : 0.0;

    AppendDouble(buffer, center.x);
    buffer.push_back(';');
    AppendDouble(buffer, center.y);
    buffer.push_back(';');
    AppendDouble(buffer, center.z);
    buffer.push_back(';');
    AppendIntegral(buffer, particle_count);
    buffer.push_back(';');
    AppendDouble(buffer, density_value);
    buffer.push_back(';');
    AppendDouble(buffer, moments.temperature);
    buffer.push_back(';');
    AppendDouble(buffer, velocity.x);
    buffer.push_back(';');
    AppendDouble(buffer, velocity.y);
    buffer.push_back(';');
    AppendDouble(buffer, velocity.z);
    buffer.push_back(';');
    AppendDouble(buffer, moments.energy);
    buffer.push_back('\n');
  }

  file.write(buffer.data(), static_cast<std::streamsize>(buffer.size()));
  file.flush();
  return file.good();
}

bool SnapshotWriter::WriteBinarySnapshot(std::ofstream& file) {
  static constexpr char kMagic[4] = {'M', 'C', '3', 'D'};
  file.write(kMagic, sizeof(kMagic));

  const std::uint32_t version = 1;
  file.write(reinterpret_cast<const char*>(&version), sizeof(version));

  const std::uint64_t record_count = static_cast<std::uint64_t>(cluster_.cells.size());
  file.write(reinterpret_cast<const char*>(&record_count),
             sizeof(record_count));

  const double density_reference = cluster_.density;
  file.write(reinterpret_cast<const char*>(&density_reference),
             sizeof(density_reference));

  if (!file) {
    return false;
  }

  for (auto& cell : cluster_.cells) {
    const auto center = cell.GetCenter();
    const auto moments = cell.ComputeMoments();
    const auto velocity = moments.velocity;
    const std::uint32_t particle_count =
        static_cast<std::uint32_t>(cell.GetParticleCount());
    const double density_value =
        cell.GetVolume() > 0 && cluster_.density != 0.0
            ? static_cast<double>(particle_count) / (cell.GetVolume() * cluster_.density)
            : 0.0;
    const double temperature = moments.temperature;
    const double energy = moments.energy;

    file.write(reinterpret_cast<const char*>(&center.x), sizeof(center.x));
    file.write(reinterpret_cast<const char*>(&center.y), sizeof(center.y));
    file.write(reinterpret_cast<const char*>(&center.z), sizeof(center.z));
    file.write(reinterpret_cast<const char*>(&particle_count),
               sizeof(particle_count));
    file.write(reinterpret_cast<const char*>(&density_value),
               sizeof(density_value));
    file.write(reinterpret_cast<const char*>(&temperature),
               sizeof(temperature));
    file.write(reinterpret_cast<const char*>(&velocity.x), sizeof(velocity.x));
    file.write(reinterpret_cast<const char*>(&velocity.y), sizeof(velocity.y));
    file.write(reinterpret_cast<const char*>(&velocity.z), sizeof(velocity.z));
    file.write(reinterpret_cast<const char*>(&energy), sizeof(energy));

    if (!file) {
      return false;
    }
  }

  file.flush();
  return file.good();
}

bool SnapshotWriter::WriteSpeedFile(const char* file_name) {
  std::ofstream file1(file_name);
  Point ap;
  for (auto& cell : cluster_.cells) {
    ap = cell.GetCenter();
    file1 << ap << cell.ComputeMoments().velocity << '\n';
  }
  file1.close();
  return file1.good();
}

bool SnapshotWriter::WriteSpeedFile() {
  std::ofstream file1("speed.dat");
  Point ap;
  for (auto& cell : cluster_.cells) {
    ap = cell.GetCenter();

    file1 << ap << "\t" << cell.ComputeMoments().velocity << '\n';
  }

  file1.close();
  return file1.good();
}

bool SnapshotWriter::WriteTimes() {
  if (cluster_.proc_id == 0) {
    std::ofstream times_file("time.dat");

    deque<double>::iterator time_iter = cluster_.times.begin();

    int i = 0;

    while (time_iter != cluster_.times.end()) {
      times_file << i << "\t" << cluster_.times[i] << "\t" << cluster_.sort_times[i] << "\t"
                 << cluster_.send_times[i] << "\t" << cluster_.bound_times[i] << "\t"
                 << cluster_.calc_times[i] << endl;

      i++;
      time_iter++;
    }

    times_file.close();
    return times_file.good();
  }
  return true;
}

bool SnapshotWriter::WriteCellFile(const std::string& init_file) {
  std::ofstream file(init_file);
  if (!file.is_open()) {
    LOG_ERROR() << "Can't open file: " << init_file;
    return false;
  }

  file << "cluster" << std::endl;

  for (auto& cell : cluster_.cells) {
    const auto moments = cell.ComputeMoments();
    Point V = moments.velocity;
    double S = std::sqrt((V.x * V.x + V.y * V.y) / 2 * moments.temperature);
    double alpha = M_PI / 2 - std::atan(V.x / V.y);
    if (V.x > 100000 * V.y) {
      V.x > 0 ? alpha = 0 : alpha = 3.141592654;
    }

    file << "\tcell" << endl;

    file << "\t\tapex\t" << cell.GetApex() << endl;

    file << "\t\tsize\t" << cell.GetSize() << endl;

    file << "\t\tS\t" << S << endl;

    file << "\t\talpha\t" << alpha << endl;

    file << "\t\tparticles\t" << cell.GetParticleCount() << endl;

    file << "\t\ttemperature\t" << moments.temperature << endl;

    file << "\t\tL\t" << cell.GetCharacteristicLength() << endl;

    file << "\tendcell" << endl;

  }

  file.flush();
  return file.good();
}

void SimulationRunner::SetSnapshotInterval(double interval) {
  if (interval > 0.0) {
    snapshot_interval_ = interval;
    next_snapshot_time_ = interval;
  } else {
    snapshot_interval_ = 0.0;
    next_snapshot_time_ = 0.0;
  }
}

bool SimulationRunner::UpdateSnapshotSchedule() {
  if (snapshot_interval_ <= 0.0) {
    return true;
  }

  if (next_snapshot_time_ <= 0.0) {
    next_snapshot_time_ = snapshot_interval_;
  }

  if (cluster_.t + kSnapshotTolerance >= next_snapshot_time_) {
    while (cluster_.t + kSnapshotTolerance >= next_snapshot_time_) {
      next_snapshot_time_ += snapshot_interval_;
    }
    return true;
  }

  return false;
}

void SimulationRunner::Compute() {
  size_t step = 0;
  while (cluster_.t < cluster_.t_end) {
    LOG_INFO() << "Time step:" << ++step;
    utils::Stopwatch stopwatch("ClusterStep");
    Step();
  }
}

void SimulationRunner::SetForceReference(double area, double speed, double alpha) {
  const double q = 0.5 * speed * speed;
  const double scale = q * area;
  if (!(area > 0) || !std::isfinite(area) || !(speed > 0) ||
      !std::isfinite(speed) || !std::isfinite(alpha) ||
      !(scale > 0) || !std::isfinite(scale))
    throw std::invalid_argument("Force coefficients require positive finite area and speed, and finite flow angle");
  if (body_loads_.is_open())
    throw std::logic_error("Cannot change force reference after recording loads");
  reference_area_ = area;
  dynamic_pressure_ = q;
  flow_angle_ = alpha;
  force_scale_ = scale;
}

void SimulationRunner::WriteBodyLoads(const StepResult& result) {
  if (!cluster_.body_ || cluster_.body_->PolygonCount() == 0) return;
  if (!body_loads_.is_open()) {
    body_loads_.open("body_loads.csv");
    body_loads_ << std::setprecision(17);
    body_loads_ << "time;dt;Fx;Fy;Fz;mean_Fx;mean_Fy;mean_Fz";
    if (force_scale_)
      body_loads_ << ";reference_area;dynamic_pressure;flow_angle;Cx;Cy;Cz;Cd;Cl;mean_Cx;mean_Cy;mean_Cz;mean_Cd;mean_Cl";
    body_loads_ << '\n';
  }
  body_loads_ << result.time << ';' << result.dt;
  for (const auto& f : {result.body_force, result.mean_body_force})
    body_loads_ << ';' << f.x << ';' << f.y << ';' << f.z;
  if (force_scale_) {
    body_loads_ << ';' << reference_area_ << ';' << dynamic_pressure_ << ';' << flow_angle_;
    for (const auto& f : {result.body_force, result.mean_body_force}) {
      const Point c = f / *force_scale_;
      const double cd = c.x * std::cos(flow_angle_) + c.y * std::sin(flow_angle_);
      const double cl = -c.x * std::sin(flow_angle_) + c.y * std::cos(flow_angle_);
      body_loads_ << ';' << c.x << ';' << c.y << ';' << c.z << ';' << cd << ';' << cl;
    }
  }
  body_loads_ << '\n';
  body_loads_.flush();
  RequireWrite(bool(body_loads_));
}

StepResult SimulationRunner::Step() {
  const auto result = cluster_.TimeStep();
  WriteBodyLoads(result);
  if (result.inner_particles > 0)
    LOG_WARNING() << "Particles inside body at time " << result.time << ": " << result.inner_particles;
  if (UpdateSnapshotSchedule()) {
    std::ostringstream name;
    name << "data" << result.time << ".dat";
    RequireWrite(writer_.WriteFile(name.str()));
  }
  return result;
}
void SimulationRunner::RequireWrite(bool success) {
  if (!success) throw std::runtime_error("Failed to write simulation output");
}
} // namespace mc3d
