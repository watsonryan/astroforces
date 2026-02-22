/**
 * @file erp_cli.cpp
 * @brief Single-state Earth radiation pressure perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/erp/erp_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 8 || argc > 11) {
    spdlog::error(
        "usage: erp_cli <x_frame_m> <y_frame_m> <z_frame_m> <vx_frame_mps> <vy_frame_mps> <vz_frame_mps> <epoch_utc_s> [mass_kg] [area_m2] [cr]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECI;

  const double mass_kg = (argc >= 9) ? std::atof(argv[8]) : 600.0;
  const double area_m2 = (argc >= 10) ? std::atof(argv[9]) : 4.0;
  const double cr = (argc >= 11) ? std::atof(argv[10]) : 1.3;

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  const astroforces::erp::ErpAccelerationModel erp{};
  const auto out = erp.evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("erp evaluation failed: status={}", static_cast<int>(out.status));
    return 2;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("p_pa={} r_earth_m={} area={} cr={}\n", out.earth_radiation_pressure_pa, out.earth_distance_m, out.area_m2, out.cr);
  return 0;
}

