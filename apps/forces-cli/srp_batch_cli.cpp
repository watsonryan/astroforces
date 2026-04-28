/**
 * @file srp_batch_cli.cpp
 * @brief Batch SRP perturbation CLI.
 * @author Watson
 */

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/srp/srp_model.hpp"

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

}  // namespace

int main(int argc, char** argv) {
  std::filesystem::path input_csv{};
  std::filesystem::path output_csv{};
  std::filesystem::path eph_file{};
  double mass_kg{600.0};
  double area_m2{4.0};
  double cr{1.3};
  int use_eclipse_flag{0};

  CLI::App app{"Batch SRP perturbation CLI"};
  app.add_option("input_csv", input_csv)->required();
  app.add_option("output_csv", output_csv)->required();
  app.add_option("jpl_ephemeris_file", eph_file)->required();
  app.add_option("mass_kg", mass_kg)->capture_default_str();
  app.add_option("area_m2", area_m2)->capture_default_str();
  app.add_option("cr", cr)->capture_default_str();
  app.add_option("use_eclipse", use_eclipse_flag)->check(CLI::Range(0, 1))->capture_default_str();
  app.footer("input row: epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps");
  CLI11_PARSE(app, argc, argv);

  const bool use_eclipse = (use_eclipse_flag != 0);

  if (!std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

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
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  auto srp = astroforces::forces::SrpAccelerationModel::Create(
      {.ephemeris_file = eph_file, .use_eclipse = use_eclipse});

  out << "epoch_utc_s,ax_mps2,ay_mps2,az_mps2,amag_mps2,solar_pressure_pa,eclipse_factor,sun_distance_m,area_m2,cr,eclipsed,status\n";

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

    const auto r = srp->evaluate(state, sc);
    out << fmt::format("{:.6f},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{},{}\n",
                       row.epoch_utc_s, r.acceleration_mps2.x, r.acceleration_mps2.y, r.acceleration_mps2.z,
                       magnitude(r.acceleration_mps2), r.solar_pressure_pa, r.eclipse_factor, r.sun_distance_m, r.area_m2, r.cr,
                       r.eclipsed ? 1 : 0, static_cast<int>(r.status));
  }

  spdlog::info("wrote srp batch output: {}", output_csv.string());
  return 0;
}
