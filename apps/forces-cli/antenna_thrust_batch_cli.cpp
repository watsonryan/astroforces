/**
 * @file antenna_thrust_batch_cli.cpp
 * @brief Batch antenna thrust perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_model.hpp"

namespace {

struct SampleRow {
  double epoch_utc_s{};
  astroforces::core::Vec3 position_eci_m{};
  astroforces::core::Vec3 velocity_eci_mps{};
};

bool parse_sample_row(const std::string& line, SampleRow& out) {
  std::stringstream ss(line);
  std::string tok;
  std::vector<double> values;
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) {
      return false;
    }
    char* end = nullptr;
    const double v = std::strtod(tok.c_str(), &end);
    if (end == tok.c_str() || *end != '\0') {
      return false;
    }
    values.push_back(v);
  }
  if (values.size() != 7U) {
    return false;
  }
  out.epoch_utc_s = values[0];
  out.position_eci_m = astroforces::core::Vec3{values[1], values[2], values[3]};
  out.velocity_eci_mps = astroforces::core::Vec3{values[4], values[5], values[6]};
  return true;
}

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

astroforces::forces::AntennaThrustDirectionMode parse_mode(const std::string& s, bool* ok) {
  if (s == "velocity") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::Velocity;
  }
  if (s == "nadir") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::Nadir;
  }
  if (s == "custom_eci") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::CustomEci;
  }
  if (s == "body_fixed") {
    *ok = true;
    return astroforces::forces::AntennaThrustDirectionMode::BodyFixed;
  }
  *ok = false;
  return astroforces::forces::AntennaThrustDirectionMode::Velocity;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3 || argc > 10) {
    spdlog::error(
        "usage: antenna_thrust_batch_cli <input_csv> <output_csv> [mass_kg] [transmit_power_w] [efficiency] [mode:velocity|nadir|custom_eci|body_fixed] [dir_x] [dir_y] [dir_z]");
    spdlog::error("input row: epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps");
    return 1;
  }

  const std::filesystem::path input_csv = argv[1];
  const std::filesystem::path output_csv = argv[2];
  const double mass_kg = (argc >= 4) ? std::atof(argv[3]) : 600.0;
  const double power_w = (argc >= 5) ? std::atof(argv[4]) : 20.0;
  const double efficiency = (argc >= 6) ? std::atof(argv[5]) : 1.0;
  const std::string mode_s = (argc >= 7) ? argv[6] : "velocity";

  bool mode_ok = false;
  const auto mode = parse_mode(mode_s, &mode_ok);
  if (!mode_ok) {
    spdlog::error("invalid mode: {}", mode_s);
    return 2;
  }

  const auto dir = astroforces::core::Vec3{
      (argc >= 8) ? std::atof(argv[7]) : 1.0,
      (argc >= 9) ? std::atof(argv[8]) : 0.0,
      (argc >= 10) ? std::atof(argv[9]) : 0.0,
  };

  std::ifstream in(input_csv);
  if (!in) {
    spdlog::error("failed to open input csv: {}", input_csv.string());
    return 3;
  }
  std::ofstream out(output_csv);
  if (!out) {
    spdlog::error("failed to open output csv: {}", output_csv.string());
    return 4;
  }

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = 4.0, .cd = 2.2, .cr = 1.3, .use_surface_model = false, .surfaces = {}};
  const astroforces::forces::AntennaThrustAccelerationModel model({
      .transmit_power_w = power_w,
      .efficiency = efficiency,
      .direction_mode = mode,
      .custom_direction_eci = dir,
      .body_axis = dir,
  });

  out << "epoch_utc_s,ax_mps2,ay_mps2,az_mps2,amag_mps2,thrust_n,effective_power_w,mass_kg,dir_x,dir_y,dir_z,status\n";

  std::string line;
  std::size_t line_no = 0;
  while (std::getline(in, line)) {
    ++line_no;
    if (line.empty()) {
      continue;
    }
    SampleRow row{};
    if (!parse_sample_row(line, row)) {
      if (line_no == 1 && line.find("epoch_utc_s") != std::string::npos) {
        continue;
      }
      spdlog::warn("skipping malformed row {}", line_no);
      continue;
    }

    astroforces::core::StateVector state{};
    state.epoch.utc_seconds = row.epoch_utc_s;
    state.position_m = row.position_eci_m;
    state.velocity_mps = row.velocity_eci_mps;
    state.frame = astroforces::core::Frame::ECI;

    const auto r = model.evaluate(state, sc);
    out << fmt::format(
        "{:.6f},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{}\n",
        row.epoch_utc_s,
        r.acceleration_mps2.x,
        r.acceleration_mps2.y,
        r.acceleration_mps2.z,
        magnitude(r.acceleration_mps2),
        r.thrust_n,
        r.effective_power_w,
        r.mass_kg,
        r.direction_eci.x,
        r.direction_eci.y,
        r.direction_eci.z,
        static_cast<int>(r.status));
  }

  spdlog::info("wrote antenna_thrust batch output: {}", output_csv.string());
  return 0;
}
