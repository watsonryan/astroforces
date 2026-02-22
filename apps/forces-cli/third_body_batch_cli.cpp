/**
 * @file third_body_batch_cli.cpp
 * @brief Batch third-body perturbation CLI.
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

#include "astroforces/forces/gravity/third_body.hpp"

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
  if (argc < 4 || argc > 6) {
    spdlog::error("usage: third_body_batch_cli <input_csv> <output_csv> <jpl_ephemeris_file> [use_sun:0|1] [use_moon:0|1]");
    spdlog::error("input row: epoch_utc_s,x_eci_m,y_eci_m,z_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps");
    return 1;
  }

  const std::filesystem::path input_csv = argv[1];
  const std::filesystem::path output_csv = argv[2];
  const std::filesystem::path eph_file = argv[3];
  const bool use_sun = (argc >= 5) ? (std::atoi(argv[4]) != 0) : true;
  const bool use_moon = (argc >= 6) ? (std::atoi(argv[5]) != 0) : true;

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

  const auto sun = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = false, .name = "third_body_sun"});
  const auto moon = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = false, .use_moon = use_moon, .name = "third_body_moon"});
  const auto total = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = use_moon, .name = "third_body_total"});

  out << "epoch_utc_s,sun_ax_mps2,sun_ay_mps2,sun_az_mps2,sun_mag_mps2,"
         "moon_ax_mps2,moon_ay_mps2,moon_az_mps2,moon_mag_mps2,"
         "total_ax_mps2,total_ay_mps2,total_az_mps2,total_mag_mps2,status\n";

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

    const auto req = astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr};
    const auto sun_c = sun->evaluate(req);
    const auto moon_c = moon->evaluate(req);
    const auto total_c = total->evaluate(req);

    astroforces::core::Status status = astroforces::core::Status::Ok;
    if (sun_c.status != astroforces::core::Status::Ok) {
      status = sun_c.status;
    }
    if (status == astroforces::core::Status::Ok && moon_c.status != astroforces::core::Status::Ok) {
      status = moon_c.status;
    }
    if (status == astroforces::core::Status::Ok && total_c.status != astroforces::core::Status::Ok) {
      status = total_c.status;
    }

    out << fmt::format(
        "{:.6f},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{}\n",
        row.epoch_utc_s, sun_c.acceleration_mps2.x, sun_c.acceleration_mps2.y, sun_c.acceleration_mps2.z,
        magnitude(sun_c.acceleration_mps2), moon_c.acceleration_mps2.x, moon_c.acceleration_mps2.y, moon_c.acceleration_mps2.z,
        magnitude(moon_c.acceleration_mps2), total_c.acceleration_mps2.x, total_c.acceleration_mps2.y, total_c.acceleration_mps2.z,
        magnitude(total_c.acceleration_mps2), static_cast<int>(status));
  }

  spdlog::info("wrote third-body batch output: {}", output_csv.string());
  return 0;
}

