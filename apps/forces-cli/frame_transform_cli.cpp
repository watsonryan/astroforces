/**
 * @file frame_transform_cli.cpp
 * @brief Frame transform utility for external cross-validation (e.g., Astropy).
 * @author Watosn
 */

#include <cstdlib>
#include <string>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/core/transforms.hpp"

namespace {

void print_vec(const char* key, const astroforces::core::Vec3& v) {
  fmt::print("{}={:.17e},{:.17e},{:.17e}\n", key, v.x, v.y, v.z);
}

}  // namespace

int main(int argc, char** argv) {
  double jd_utc{};
  double jd_tt{};
  double cip_x_rad{};
  double cip_y_rad{};
  double cip_s_rad{};
  double xp_rad{};
  double yp_rad{};
  double dut1_s{};
  double lod_s{};
  double dX_rad{};
  double dY_rad{};
  double rx_gcrf_m{};
  double ry_gcrf_m{};
  double rz_gcrf_m{};
  double vx_gcrf_mps{};
  double vy_gcrf_mps{};
  double vz_gcrf_mps{};
  double cip_x_dot_rad_s{};
  double cip_y_dot_rad_s{};
  double cip_s_dot_rad_s{};
  double xp_dot_rad_s{};
  double yp_dot_rad_s{};
  double dX_dot_rad_s{};
  double dY_dot_rad_s{};

  CLI::App app{"Frame transform utility for external cross-validation"};
  app.require_subcommand(1);

  auto* gcrf_to_itrf = app.add_subcommand("gcrf_to_itrf", "Transform GCRF position/velocity to ITRF and back");
  gcrf_to_itrf->add_option("jd_utc", jd_utc)->required();
  gcrf_to_itrf->add_option("jd_tt", jd_tt)->required();
  gcrf_to_itrf->add_option("cip_x_rad", cip_x_rad)->required();
  gcrf_to_itrf->add_option("cip_y_rad", cip_y_rad)->required();
  gcrf_to_itrf->add_option("cip_s_rad", cip_s_rad)->required();
  gcrf_to_itrf->add_option("xp_rad", xp_rad)->required();
  gcrf_to_itrf->add_option("yp_rad", yp_rad)->required();
  gcrf_to_itrf->add_option("dut1_s", dut1_s)->required();
  gcrf_to_itrf->add_option("lod_s", lod_s)->required();
  gcrf_to_itrf->add_option("dX_rad", dX_rad)->required();
  gcrf_to_itrf->add_option("dY_rad", dY_rad)->required();
  gcrf_to_itrf->add_option("rx_gcrf_m", rx_gcrf_m)->required();
  gcrf_to_itrf->add_option("ry_gcrf_m", ry_gcrf_m)->required();
  gcrf_to_itrf->add_option("rz_gcrf_m", rz_gcrf_m)->required();
  gcrf_to_itrf->add_option("vx_gcrf_mps", vx_gcrf_mps)->required();
  gcrf_to_itrf->add_option("vy_gcrf_mps", vy_gcrf_mps)->required();
  gcrf_to_itrf->add_option("vz_gcrf_mps", vz_gcrf_mps)->required();
  gcrf_to_itrf->add_option("--cip-x-dot-rad-s", cip_x_dot_rad_s);
  gcrf_to_itrf->add_option("--cip-y-dot-rad-s", cip_y_dot_rad_s);
  gcrf_to_itrf->add_option("--cip-s-dot-rad-s", cip_s_dot_rad_s);
  gcrf_to_itrf->add_option("--xp-dot-rad-s", xp_dot_rad_s);
  gcrf_to_itrf->add_option("--yp-dot-rad-s", yp_dot_rad_s);
  gcrf_to_itrf->add_option("--dX-dot-rad-s", dX_dot_rad_s);
  gcrf_to_itrf->add_option("--dY-dot-rad-s", dY_dot_rad_s);
  CLI11_PARSE(app, argc, argv);

  astroforces::core::CelestialIntermediatePole cip{};
  cip.x_rad = cip_x_rad;
  cip.y_rad = cip_y_rad;
  cip.s_rad = cip_s_rad;

  astroforces::core::EarthOrientation eop{};
  eop.xp_rad = xp_rad;
  eop.yp_rad = yp_rad;
  eop.dut1_s = dut1_s;
  eop.lod_s = lod_s;
  eop.dX_rad = dX_rad;
  eop.dY_rad = dY_rad;
  astroforces::core::CelestialIntermediatePoleRate cip_rate{};
  astroforces::core::EarthOrientationRate eop_rate{};
  const bool have_rates =
      gcrf_to_itrf->count("--cip-x-dot-rad-s") || gcrf_to_itrf->count("--cip-y-dot-rad-s")
      || gcrf_to_itrf->count("--cip-s-dot-rad-s") || gcrf_to_itrf->count("--xp-dot-rad-s")
      || gcrf_to_itrf->count("--yp-dot-rad-s") || gcrf_to_itrf->count("--dX-dot-rad-s")
      || gcrf_to_itrf->count("--dY-dot-rad-s");
  if (have_rates) {
    cip_rate.x_rad_s = cip_x_dot_rad_s;
    cip_rate.y_rad_s = cip_y_dot_rad_s;
    cip_rate.s_rad_s = cip_s_dot_rad_s;
    eop_rate.xp_rad_s = xp_dot_rad_s;
    eop_rate.yp_rad_s = yp_dot_rad_s;
    eop_rate.dX_rad_s = dX_dot_rad_s;
    eop_rate.dY_rad_s = dY_dot_rad_s;
  }

  const astroforces::core::Vec3 r_gcrf{rx_gcrf_m, ry_gcrf_m, rz_gcrf_m};
  const astroforces::core::Vec3 v_gcrf{vx_gcrf_mps, vy_gcrf_mps, vz_gcrf_mps};

  astroforces::core::RotationWithDerivative rd{};
  if (have_rates) {
    rd = astroforces::core::gcrf_to_itrf_rotation_with_derivative_exact(jd_utc, jd_tt, cip, cip_rate, eop, eop_rate);
  } else {
    rd = astroforces::core::gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  }
  const auto r_itrf = astroforces::core::mat_vec(rd.r, r_gcrf);
  const auto v_itrf = astroforces::core::mat_vec(rd.r, v_gcrf) + astroforces::core::mat_vec(rd.dr, r_gcrf);
  const auto rt = astroforces::core::mat_transpose(rd.r);
  const auto r_gcrf_back = astroforces::core::mat_vec(rt, r_itrf);
  const auto v_gcrf_back = astroforces::core::mat_vec(rt, v_itrf - astroforces::core::mat_vec(rd.dr, r_gcrf_back));

  print_vec("r_itrf_m", r_itrf);
  print_vec("v_itrf_mps", v_itrf);
  print_vec("r_gcrf_back_m", r_gcrf_back);
  print_vec("v_gcrf_back_mps", v_gcrf_back);
  return 0;
}
