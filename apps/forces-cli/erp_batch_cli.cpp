/**
 * @file erp_batch_cli.cpp
 * @brief Batch Earth radiation pressure perturbation CLI.
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

#include "astroforces/forces/surface/erp/erp_model.hpp"

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
  if (argc < 3 || argc > 6) {
    spdlog::error("usage: erp_batch_cli <input_csv> <output_csv> [mass_kg] [area_m2] [cr]");
    spdlog::error("input row: epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps");
    return 1;
  }

  const std::filesystem::path input_csv = argv[1];
  const std::filesystem::path output_csv = argv[2];
  const double mass_kg = (argc >= 4) ? std::atof(argv[3]) : 600.0;
  const double area_m2 = (argc >= 5) ? std::atof(argv[4]) : 4.0;
  const double cr = (argc >= 6) ? std::atof(argv[5]) : 1.3;

  std::ifstream in(input_csv);
  if (!in) {
    spdlog::error("failed to open input csv: {}", input_csv.string());
    return 2;
  }
  std::ofstream out(output_csv);
  if (!out) {
    spdlog::error("failed to open output csv: {}", output_csv.string());
    return 3;
  }

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};
  const astroforces::erp::ErpAccelerationModel erp{};

  out << "epoch_utc_s,ax_mps2,ay_mps2,az_mps2,amag_mps2,earth_radiation_pressure_pa,earth_distance_m,area_m2,cr,status\n";

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

    const auto r = erp.evaluate(state, sc);
    out << fmt::format("{:.6f},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{}\n",
                       row.epoch_utc_s, r.acceleration_mps2.x, r.acceleration_mps2.y, r.acceleration_mps2.z,
                       magnitude(r.acceleration_mps2), r.earth_radiation_pressure_pa, r.earth_distance_m, r.area_m2, r.cr,
                       static_cast<int>(r.status));
  }

  spdlog::info("wrote erp batch output: {}", output_csv.string());
  return 0;
}
