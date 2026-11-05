#include "statistical_diagnostics.h"

#include <filesystem>
#include <iomanip>
#include <sstream>

#include "../third_party/json.hpp"
namespace mc3d {
namespace {
void AtomicWrite(const std::string& path, const std::string& text) {
  std::ofstream out(path + ".tmp");
  out << text;
  out.close();
  if (!out) throw std::runtime_error("Cannot write statistics: " + path);
  std::filesystem::rename(path + ".tmp", path);
}
void Vector(std::ostream& out, Point p) {
  out << ';' << p.x << ';' << p.y << ';' << p.z;
}
}  // namespace
ConservedQuantity StatisticalDiagnostics::Inventory() const {
  ConservedQuantity result;
  for (const auto& cell : cluster_.cells)
    for (const auto& p : cell.particles) result.Add(p.velocity);
  return result;
}
double StatisticalDiagnostics::WallEnergy() const {
  double result = 0;
  if (cluster_.body_ && cluster_.body_->SurfaceHeatingEnabled())
    for (size_t i = 0; i < cluster_.body_->PolygonCount(); ++i) {
      const auto& p = cluster_.body_->GetPolygon(i);
      result += p.heat_capacity * p.wall_temperature;
    }
  return result;
}
StatisticalDiagnostics::StatisticalDiagnostics(CellCluster& cluster,
                                               const SimulationConfig& config)
    : cluster_(cluster),
      statistics_(config.statistics_enabled),
      balance_(config.balance_enabled),
      weight_(1 / cluster.density) {
  if (!(weight_ > 0) || !std::isfinite(weight_))
    throw std::invalid_argument(
        "Statistics require positive reference particle density");
  fixed_temperature_ = config.wall_temperature;
  auto intervals = config.statistics_windows;
  if (intervals.empty()) intervals.push_back({0, config.end_time});
  if (statistics_)
    for (auto w : intervals)
      windows_.push_back(
          {w[0], w[1], 0, std::vector<WindowMoments>(cluster.cells.size()),
           std::vector<SurfaceSum>(cluster.body_ ? cluster.body_->PolygonCount()
                                                 : 0)});
  initial_ = last_ = Inventory();
  wall_initial_ = WallEnergy();
  if (balance_)
    AtomicWrite(
        "conservation.csv",
        "step;time;end_time;dt;N_before;N_after;N_in;N_out;M_before;M_after;M_"
        "in;M_out;Px_before;Py_before;Pz_before;Px_after;Py_after;Pz_after;Px_"
        "in;Py_in;Pz_in;Px_out;Py_out;Pz_out;E_before;E_after;E_in;E_out;Jbody_"
        "x;Jbody_y;Jbody_z;Jmirror_x;Jmirror_y;Jmirror_z;Qbody;Qmirror;delta_"
        "Uwall;Qthermostat;Rcount;Rmass;Rp_x;Rp_y;Rp_z;RE;Rwall\n");
}
void StatisticalDiagnostics::BeforeStep() {
  start_ = cluster_.t;
  dt_ = cluster_.dt;
  if (balance_) {
    before_ = Inventory();
    wall_before_ = WallEnergy();
  }
  if (!statistics_) return;
  for (auto& w : windows_) {
    const double h =
        std::max(0., std::min(start_ + dt_, w.end) - std::max(start_, w.start));
    if (h == 0) continue;
    w.elapsed += h;
    for (size_t i = 0; i < cluster_.cells.size(); ++i) {
      auto& sum = w.cells[i];
      sum.duration += h;
      for (const auto& p : cluster_.cells[i].particles) sum.Add(p.velocity, h);
    }
  }
}
void StatisticalDiagnostics::AfterStep() {
  ++steps_;
  if (statistics_)
    for (auto& w : windows_) {
      const double h = std::max(
          0., std::min(start_ + dt_, w.end) - std::max(start_, w.start));
      if (h == 0) continue;
      for (size_t i = 0; i < w.surface.size(); ++i) {
        const auto& p = cluster_.body_->GetPolygon(i);
        w.surface[i].impulse += p.step_impulse * (weight_ * h / dt_);
        w.surface[i].heat += p.step_heat * (weight_ * h / dt_);
      }
    }
  if (balance_) {
    last_ = Inventory();
    const auto& in = cluster_.step_incoming;
    const auto& out = cluster_.step_outgoing;
    const auto jp = cluster_.step_body_impulse * weight_,
               jm = cluster_.step_mirror_impulse * weight_;
    const double qb = cluster_.step_body_heat * weight_,
                 qm = cluster_.step_mirror_heat * weight_;
    const auto rn = last_.count - before_.count - in.count + out.count;
    const auto rp =
        (last_.momentum - before_.momentum - in.momentum + out.momentum) *
            weight_ +
        jp + jm;
    const double re =
        (last_.energy - before_.energy - in.energy + out.energy) * weight_ +
        qb + qm;
    const bool dynamic =
        cluster_.body_ && cluster_.body_->SurfaceHeatingEnabled();
    const double du = WallEnergy() - wall_before_;
    if (dynamic)
      max_wall_residual_ = std::max(max_wall_residual_, std::abs(du - qb));
    incoming_ += in;
    outgoing_ += out;
    body_impulse_ += jp;
    mirror_impulse_ += jm;
    body_heat_ += qb;
    mirror_heat_ += qm;
    count_residual_ += rn;
    momentum_residual_ += rp;
    energy_residual_ += re;
    max_energy_residual_ = std::max(max_energy_residual_, std::abs(re));
    max_momentum_residual_ = std::max(max_momentum_residual_, rp.Mod());
    std::ofstream f("conservation.csv", std::ios::app);
    f << std::setprecision(17);
    f << steps_ << ';' << start_ << ';' << cluster_.t << ';' << dt_ << ';'
      << before_.count << ';' << last_.count << ';' << in.count << ';'
      << out.count;
    for (auto n : {before_.count, last_.count, in.count, out.count})
      f << ';' << n * weight_;
    for (auto p : {before_.momentum, last_.momentum, in.momentum, out.momentum})
      Vector(f, p * weight_);
    for (auto e : {before_.energy, last_.energy, in.energy, out.energy})
      f << ';' << e * weight_;
    Vector(f, jp);
    Vector(f, jm);
    f << ';' << qb << ';' << qm << ';' << du << ';' << (dynamic ? 0 : qb) << ';'
      << rn << ';' << rn * weight_;
    Vector(f, rp);
    f << ';' << re << ';' << (dynamic ? du - qb : 0) << '\n';
    f.close();
    if (!f) throw std::runtime_error("Cannot write conservation.csv");
  }
  for (const auto& w : windows_)
    if (start_ < w.end && cluster_.t >= w.end) {
      Write();
      break;
    }
}
void StatisticalDiagnostics::Write() {
  if (statistics_)
    for (size_t k = 0; k < windows_.size(); ++k) {
      const auto& w = windows_[k];
      std::ostringstream f;
      f << std::setprecision(17);
      f << "window;start;end;observed_time;complete;cell;status;Vgas;A0;origin_"
           "x;origin_y;origin_z;A1x;A1y;A1z;A2xx;A2xy;A2xz;A2yy;A2yz;A2zz;A3x;"
           "A3y;A3z;mean_N;density;ux;uy;uz;Tstat;Pxx;Pxy;Pxz;Pyy;Pyz;Pzz;qx;"
           "qy;qz\n";
      for (size_t i = 0; i < w.cells.size(); ++i) {
        const auto& a = w.cells[i];
        const auto& cell = cluster_.cells[i];
        const double volume = cell.GetVolume();
        const char* status = volume == 0      ? "solid"
                             : w.elapsed == 0 ? "unobserved"
                             : a.a0 == 0      ? "empty"
                                              : "valid";
        f << k << ';' << w.start << ';' << w.end << ';' << w.elapsed << ';'
          << (cluster_.t >= w.end) << ';' << i << ';' << status << ';' << volume
          << ';' << a.a0;
        for (auto x : a.origin) f << ';' << x;
        for (auto x : a.a1) f << ';' << x;
        for (int j : {0, 1, 2, 4, 5, 8}) f << ';' << a.a2[j];
        for (auto x : a.a3) f << ';' << x;
        if (w.elapsed > 0 && volume > 0)
          f << ';' << a.a0 / w.elapsed << ';'
            << weight_ * a.a0 / (volume * w.elapsed);
        else
          f << ";;";
        if (a.a0 > 0 && volume > 0) {
          const auto r = a.Evaluate(weight_, volume);
          Vector(f, r.velocity);
          f << ';' << r.temperature;
          for (int j : {0, 1, 2, 4, 5, 8}) f << ';' << r.pressure[j];
          Vector(f, r.heat);
        } else
          f << ";;;;;;;;;;;;;";
        f << '\n';
      }
      AtomicWrite("cell_statistics_" + std::to_string(k) + ".csv", f.str());
      std::ostringstream s;
      s << std::setprecision(17)
        << "window;start;end;observed_time;complete;polygon;area;wall_"
           "temperature_at_output;heat_capacity;Jx;Jy;Jz;Q;Fx;Fy;Fz;power;heat_"
           "flux\n";
      for (size_t i = 0; i < w.surface.size(); ++i) {
        const auto& p = cluster_.body_->GetPolygon(i);
        const auto& a = w.surface[i];
        const double area = .5 * (p.p2 - p.p1).Cross(p.p3 - p.p1).Mod();
        s << k << ';' << w.start << ';' << w.end << ';' << w.elapsed << ';'
          << (cluster_.t >= w.end) << ';' << i << ';' << area << ';'
          << (cluster_.body_->SurfaceHeatingEnabled() ? p.wall_temperature
                                                      : fixed_temperature_)
          << ';' << p.heat_capacity;
        Vector(s, a.impulse);
        s << ';' << a.heat;
        if (w.elapsed > 0) {
          Vector(s, a.impulse / w.elapsed);
          s << ';' << a.heat / w.elapsed << ';' << a.heat / (w.elapsed * area);
        } else
          s << ";;;;;";
        s << '\n';
      }
      AtomicWrite("surface_statistics_" + std::to_string(k) + ".csv", s.str());
    }
  if (balance_) {
    auto quantity = [&](const ConservedQuantity& x) {
      return nlohmann::json{
          {"N", x.count},
          {"mass", x.count * weight_},
          {"momentum",
           {x.momentum.x * weight_, x.momentum.y * weight_,
            x.momentum.z * weight_}},
          {"energy", x.energy * weight_},
          {"absolute_momentum", x.absolute_momentum * weight_}};
    };
    nlohmann::json j = {
        {"schema_version", 1},
        {"steps", steps_},
        {"time", cluster_.t},
        {"initial", quantity(initial_)},
        {"final", quantity(last_)},
        {"incoming", quantity(incoming_)},
        {"outgoing", quantity(outgoing_)},
        {"body_heat", body_heat_},
        {"mirror_heat", mirror_heat_},
        {"body_impulse", {body_impulse_.x, body_impulse_.y, body_impulse_.z}},
        {"mirror_impulse",
         {mirror_impulse_.x, mirror_impulse_.y, mirror_impulse_.z}},
        {"residual_count", count_residual_},
        {"residual_mass", count_residual_ * weight_},
        {"residual_momentum",
         {momentum_residual_.x, momentum_residual_.y, momentum_residual_.z}},
        {"residual_energy", energy_residual_},
        {"max_energy_residual", max_energy_residual_},
        {"max_momentum_residual", max_momentum_residual_},
        {"wall_energy_change", WallEnergy() - wall_initial_},
        {"thermostat_heat",
         cluster_.body_ && cluster_.body_->SurfaceHeatingEnabled()
             ? 0
             : body_heat_},
        {"max_wall_residual", max_wall_residual_}};
    AtomicWrite("conservation_summary.json", j.dump(2) + "\n");
  }
}
}  // namespace mc3d
