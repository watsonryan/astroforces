/**
 * @file test_gravity_sph.cpp
 * @brief Full SPH gravity model tests.
 * @author Watosn
 */

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/gravity_sph_model.hpp"

namespace {

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

std::filesystem::path write_test_gfc() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_gm.gfc";
  std::ofstream out(p);
  out << "modelname TEST\n";
  out << "earth_gravity_constant 3.986004418e14\n";
  out << "radius 6378136.3\n";
  out << "max_degree 4\n";
  out << "tide_system tide_free\n";
  out << "end_of_head\n";
  out << "gfc 0 0 1.0 0.0 0.0 0.0\n";
  out << "gfc 2 0 -4.84165143790815e-4 0.0 0.0 0.0\n";
  out << "gfc 2 1 1.0e-10 -2.0e-10 0.0 0.0\n";
  out << "gfc 2 2 2.0e-6 -1.0e-6 0.0 0.0\n";
  out << "gfc 3 0 9.57161e-7 0.0 0.0 0.0\n";
  out << "gfc 3 1 2.0e-7 1.0e-7 0.0 0.0\n";
  out << "gfc 3 2 1.0e-7 -2.0e-7 0.0 0.0\n";
  out << "gfc 3 3 5.0e-8 4.0e-8 0.0 0.0\n";
  out.close();
  return p;
}

std::filesystem::path write_test_eop() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_eop.txt";
  std::ofstream out(p);
  out << " 24 12 31 60676.00 I  0.123456 0.000001  0.234567 0.000001 I  0.1000000 0.0001000  0.900 0.010\n";
  out << " 25  1  1 60677.00 I  0.223456 0.000001  0.334567 0.000001 I  0.1200000 0.0001000  1.100 0.010\n";
  out << " 25  1  2 60678.00 P  0.323456 0.000001  0.434567 0.000001 P  0.1500000 0.0001000  1.300 0.010\n";
  out.close();
  return p;
}

std::filesystem::path write_test_constituent_tide() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_tide.txt";
  std::ofstream out(p);
  out << "header1\n";
  out << "header2\n";
  out << "header3\n";
  out << "header4\n";
  out << "055.565 TESTWAVE 2 1  1.0  0.5  0.2  -0.1\n";
  out.close();
  return p;
}

std::filesystem::path write_test_ocean_pole_tide() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_ocean_pole_tide.txt";
  std::ofstream out(p);
  out << "# n m cnmp cnmm snmp snmm\n";
  out << "2 1 1.0e-2 5.0e-3 -3.0e-3 2.0e-3\n";
  out.close();
  return p;
}

std::filesystem::path write_test_aod() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_aod.txt";
  std::ofstream out(p);
  out << "DATA SET 1: 2 COEFFICIENTS FOR 2001-09-09 01:46:39 OF TYPE glo\n";
  out << "2 0 1.0e-11 0.0\n";
  out << "2 1 2.0e-11 -1.0e-11\n";
  out << "DATA SET 2: 2 COEFFICIENTS FOR 2001-09-09 07:46:39 OF TYPE glo\n";
  out << "2 0 3.0e-11 0.0\n";
  out << "2 1 4.0e-11 -2.0e-11\n";
  out.close();
  return p;
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  const auto gravity_file = write_test_gfc();

  fs::path eph_path;
  if (const char* env = std::getenv("ASTROFORCES_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(ASTROFORCES_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
  }

  astroforces::core::StateVector state{};
  state.frame = astroforces::core::Frame::ECEF;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::core::Vec3{6878137.0, 1200.0, -3400.0};
  state.velocity_mps = astroforces::core::Vec3{0.0, 7600.0, 0.0};

  const auto central_only = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file, .max_degree = 4, .use_central = true, .use_sph = false, .use_solid_earth_tides = false});
  const auto full_sph = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file, .max_degree = 4, .use_central = true, .use_sph = true, .use_solid_earth_tides = false});

  const auto out_c = central_only->evaluate(state);
  const auto out_s = full_sph->evaluate(state);
  if (out_c.status != astroforces::core::Status::Ok || out_s.status != astroforces::core::Status::Ok) {
    spdlog::error("gravity evaluation failed");
    return 1;
  }
  if (!finite_vec(out_c.acceleration_mps2) || !finite_vec(out_s.acceleration_mps2)) {
    spdlog::error("gravity acceleration not finite");
    return 2;
  }

  const auto sph_delta = astroforces::core::Vec3{
      out_s.acceleration_mps2.x - out_c.acceleration_mps2.x,
      out_s.acceleration_mps2.y - out_c.acceleration_mps2.y,
      out_s.acceleration_mps2.z - out_c.acceleration_mps2.z};
  if (!(astroforces::core::norm(sph_delta) > 0.0)) {
    spdlog::error("sph harmonics had no effect");
    return 3;
  }

  state.frame = astroforces::core::Frame::ECI;
  const auto out_eci = full_sph->evaluate(state);
  if (out_eci.status != astroforces::core::Status::Ok || !finite_vec(out_eci.acceleration_mps2)) {
    spdlog::error("eci evaluation failed");
    return 4;
  }

  if (fs::exists(eph_path)) {
    const auto no_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .ephemeris_file = eph_path,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false});
    const auto with_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .ephemeris_file = eph_path,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = true});
    const auto out_no_tides = no_tides->evaluate(state);
    const auto out_with_tides = with_tides->evaluate(state);
    if (out_no_tides.status != astroforces::core::Status::Ok || out_with_tides.status != astroforces::core::Status::Ok) {
      spdlog::error("tide model evaluation failed");
      return 5;
    }
    const auto tide_delta = astroforces::core::Vec3{
        out_with_tides.acceleration_mps2.x - out_no_tides.acceleration_mps2.x,
        out_with_tides.acceleration_mps2.y - out_no_tides.acceleration_mps2.y,
        out_with_tides.acceleration_mps2.z - out_no_tides.acceleration_mps2.z};
    if (!(astroforces::core::norm(tide_delta) > 0.0)) {
      spdlog::error("solid Earth tides had no effect");
      return 6;
    }
  }

  const auto eop_file = write_test_eop();
  if (fs::exists(eop_file)) {
    const auto no_pole_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_pole_tide_solid = false,
         .use_pole_tide_ocean = false});
    const auto with_pole_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .eop_finals_file = eop_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_pole_tide_solid = true,
         .use_pole_tide_ocean = true});

    const auto out_no_pole = no_pole_tides->evaluate(state);
    const auto out_with_pole = with_pole_tides->evaluate(state);
    if (out_no_pole.status != astroforces::core::Status::Ok || out_with_pole.status != astroforces::core::Status::Ok) {
      spdlog::error("pole tide model evaluation failed");
      return 8;
    }

    const auto pole_delta = astroforces::core::Vec3{
        out_with_pole.acceleration_mps2.x - out_no_pole.acceleration_mps2.x,
        out_with_pole.acceleration_mps2.y - out_no_pole.acceleration_mps2.y,
        out_with_pole.acceleration_mps2.z - out_no_pole.acceleration_mps2.z};
    if (!(astroforces::core::norm(pole_delta) > 0.0)
        || !(astroforces::core::norm(out_with_pole.pole_tide_solid_mps2) > 0.0)
        || !(astroforces::core::norm(out_with_pole.pole_tide_ocean_mps2) > 0.0)) {
      spdlog::error("pole tides had no effect");
      return 9;
    }

    const auto ocean_pole_file = write_test_ocean_pole_tide();
    const auto with_pole_file = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .eop_finals_file = eop_file,
         .ocean_pole_tide_file = ocean_pole_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_pole_tide_solid = false,
         .use_pole_tide_ocean = true});
    const auto out_with_pole_file = with_pole_file->evaluate(state);
    if (out_with_pole_file.status != astroforces::core::Status::Ok) {
      spdlog::error("file-driven pole tide model evaluation failed");
      return 14;
    }
    if (!(astroforces::core::norm(out_with_pole_file.pole_tide_ocean_mps2) > 0.0)) {
      spdlog::error("file-driven ocean pole tide had no effect");
      return 15;
    }
    const auto ocean_pole_delta = astroforces::core::Vec3{
        out_with_pole_file.pole_tide_ocean_mps2.x - out_with_pole.pole_tide_ocean_mps2.x,
        out_with_pole_file.pole_tide_ocean_mps2.y - out_with_pole.pole_tide_ocean_mps2.y,
        out_with_pole_file.pole_tide_ocean_mps2.z - out_with_pole.pole_tide_ocean_mps2.z};
    if (!(astroforces::core::norm(ocean_pole_delta) > 0.0)) {
      spdlog::error("file-driven ocean pole tide did not change fallback result");
      return 16;
    }
  }

  {
    const auto no_tide2 = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = true,
         .use_sun_tide = false,
         .use_moon_tide = false,
         .use_solid_earth_tide2 = false});

    const auto with_tide2 = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = true,
         .use_sun_tide = false,
         .use_moon_tide = false,
         .use_solid_earth_tide2 = true});

    const auto out_no_tide2 = no_tide2->evaluate(state);
    const auto out_with_tide2 = with_tide2->evaluate(state);
    if (out_no_tide2.status != astroforces::core::Status::Ok || out_with_tide2.status != astroforces::core::Status::Ok) {
      spdlog::error("solid Earth tide2 evaluation failed");
      return 12;
    }
    const auto tide2_delta = astroforces::core::Vec3{
        out_with_tide2.acceleration_mps2.x - out_no_tide2.acceleration_mps2.x,
        out_with_tide2.acceleration_mps2.y - out_no_tide2.acceleration_mps2.y,
        out_with_tide2.acceleration_mps2.z - out_no_tide2.acceleration_mps2.z};
    if (!(astroforces::core::norm(tide2_delta) > 0.0)
        || !(astroforces::core::norm(out_with_tide2.solid_tide_freqdep_mps2) > 0.0)) {
      spdlog::error("solid Earth tide2 had no effect");
      return 13;
    }
  }

  const auto constituent_tide_file = write_test_constituent_tide();
  if (fs::exists(constituent_tide_file)) {
    const auto no_constituent = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_ocean_tide = false,
         .use_atmos_tide = false});

    const auto with_constituent = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .ocean_tide_file = constituent_tide_file,
         .atmos_tide_file = constituent_tide_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_ocean_tide = true,
         .use_atmos_tide = true});

    const auto out_no_constituent = no_constituent->evaluate(state);
    const auto out_with_constituent = with_constituent->evaluate(state);
    if (out_no_constituent.status != astroforces::core::Status::Ok
        || out_with_constituent.status != astroforces::core::Status::Ok) {
      spdlog::error("constituent tide model evaluation failed");
      return 10;
    }

    const auto constituent_delta = astroforces::core::Vec3{
        out_with_constituent.acceleration_mps2.x - out_no_constituent.acceleration_mps2.x,
        out_with_constituent.acceleration_mps2.y - out_no_constituent.acceleration_mps2.y,
        out_with_constituent.acceleration_mps2.z - out_no_constituent.acceleration_mps2.z};
    if (!(astroforces::core::norm(constituent_delta) > 0.0)
        || !(astroforces::core::norm(out_with_constituent.ocean_tide_mps2) > 0.0)
        || !(astroforces::core::norm(out_with_constituent.atmos_tide_mps2) > 0.0)) {
      spdlog::error("constituent tides had no effect");
      return 11;
    }
  }

  {
    const auto aod_file = write_test_aod();
    const auto no_aod = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_aod = false});

    const auto with_aod = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .aod_file = aod_file,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false,
         .use_aod = true});

    const auto out_no_aod = no_aod->evaluate(state);
    const auto out_with_aod = with_aod->evaluate(state);
    if (out_no_aod.status != astroforces::core::Status::Ok || out_with_aod.status != astroforces::core::Status::Ok) {
      spdlog::error("AOD model evaluation failed");
      return 17;
    }
    const auto aod_delta = astroforces::core::Vec3{
        out_with_aod.acceleration_mps2.x - out_no_aod.acceleration_mps2.x,
        out_with_aod.acceleration_mps2.y - out_no_aod.acceleration_mps2.y,
        out_with_aod.acceleration_mps2.z - out_no_aod.acceleration_mps2.z};
    if (!(astroforces::core::norm(aod_delta) > 0.0) || !(astroforces::core::norm(out_with_aod.aod_mps2) > 0.0)) {
      spdlog::error("AOD had no effect");
      return 18;
    }
  }

  state.frame = astroforces::core::Frame::NED;
  const auto out_bad = full_sph->evaluate(state);
  if (out_bad.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("expected invalid frame status");
    return 7;
  }

  return 0;
}
