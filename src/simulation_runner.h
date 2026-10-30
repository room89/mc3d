#pragma once
#include "cell_cluster.h"
namespace mc3d {
class SnapshotWriter {
 public:
  explicit SnapshotWriter(CellCluster& cluster) : cluster_(cluster) {}
  void SetBinaryOutput(bool enabled);
  bool WriteFile(const std::string& name);
  bool WriteFile();
  bool WriteSpeedFile(const char* name);
  bool WriteSpeedFile();
  bool WriteTimes();
  bool WriteCellFile(const std::string& name);
 private:
  CellCluster& cluster_;
  bool binary_output_ = false;
  bool WriteTextSnapshot(std::ofstream& file);
  bool WriteBinarySnapshot(std::ofstream& file);
};
class SimulationRunner {
 public:
  explicit SimulationRunner(CellCluster& cluster) : cluster_(cluster), writer_(cluster) {}
  SnapshotWriter& Writer() { return writer_; }
  void SetSnapshotInterval(double interval);
  double SnapshotInterval() const { return snapshot_interval_; }
  double NextSnapshotTime() const { return next_snapshot_time_; }
  // rho_infinity = 1 in normalized solver units.
  void SetForceReference(double area, double speed, double alpha);
  StepResult Step();
  void Compute();
  static void RequireWrite(bool success);
 private:
  CellCluster& cluster_;
  SnapshotWriter writer_;
  double snapshot_interval_ = 0;
  double next_snapshot_time_ = 0;
  bool UpdateSnapshotSchedule();
  void WriteBodyLoads(const StepResult& result);
  std::ofstream body_loads_;
  std::optional<double> force_scale_;
  double reference_area_ = 0;
  double dynamic_pressure_ = 0;
  double flow_angle_ = 0;
};
}
