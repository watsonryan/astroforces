/**
 * @file test_orbital_elements.cpp
 * @brief Classical orbital-element regression tests.
 * @author Watson
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/orbital_elements.hpp"

namespace {

using astroforces::core::Vec3;

bool approx_abs(const double a, const double b, const double tol) {
  return std::abs(a - b) <= tol;
}

double wrap_diff(double a, double b) {
  double d = std::fmod(a - b, astroforces::core::constants::kTwoPi);
  if (d > astroforces::core::constants::kPi) {
    d -= astroforces::core::constants::kTwoPi;
  } else if (d < -astroforces::core::constants::kPi) {
    d += astroforces::core::constants::kTwoPi;
  }
  return d;
}

bool approx_angle(const double a, const double b, const double tol) {
  return std::abs(wrap_diff(a, b)) <= tol;
}

struct StatePair {
  Vec3 r_m{};
  Vec3 v_mps{};
};

StatePair coe_to_state(const double a_m,
                       const double e,
                       const double i_rad,
                       const double raan_rad,
                       const double argp_rad,
                       const double nu_rad,
                       const double mu_m3_s2) {
  const double p = a_m * (1.0 - e * e);
  const double r = p / (1.0 + e * std::cos(nu_rad));
  const double h = std::sqrt(mu_m3_s2 * p);

  const Vec3 r_pf{r * std::cos(nu_rad), r * std::sin(nu_rad), 0.0};
  const Vec3 v_pf{-mu_m3_s2 / h * std::sin(nu_rad), mu_m3_s2 / h * (e + std::cos(nu_rad)), 0.0};

  const double cO = std::cos(raan_rad);
  const double sO = std::sin(raan_rad);
  const double ci = std::cos(i_rad);
  const double si = std::sin(i_rad);
  const double cw = std::cos(argp_rad);
  const double sw = std::sin(argp_rad);

  const double q11 = cO * cw - sO * sw * ci;
  const double q12 = -cO * sw - sO * cw * ci;
  const double q21 = sO * cw + cO * sw * ci;
  const double q22 = -sO * sw + cO * cw * ci;
  const double q31 = sw * si;
  const double q32 = cw * si;

  return StatePair{
      .r_m = Vec3{
          q11 * r_pf.x + q12 * r_pf.y,
          q21 * r_pf.x + q22 * r_pf.y,
          q31 * r_pf.x + q32 * r_pf.y,
      },
      .v_mps = Vec3{
          q11 * v_pf.x + q12 * v_pf.y,
          q21 * v_pf.x + q22 * v_pf.y,
          q31 * v_pf.x + q32 * v_pf.y,
      },
  };
}

}  // namespace

int main() {
  using namespace astroforces::core;

  const double mu = constants::kEarthMuM3S2;

  {
    const double radius_m = 7000e3;
    const double speed_mps = std::sqrt(mu / radius_m);
    const auto coe =
        state_eci_to_classical_orbital_elements(Vec3{radius_m, 0.0, 0.0}, Vec3{0.0, speed_mps, 0.0}, mu);
    if (coe.status != Status::Ok) {
      spdlog::error("circular-equatorial conversion failed");
      return 1;
    }
    if (!approx_abs(coe.semi_major_axis_m, radius_m, 1e-6) || !approx_abs(coe.eccentricity, 0.0, 1e-12)
        || !approx_abs(coe.inclination_rad, 0.0, 1e-12) || !approx_abs(coe.raan_rad, 0.0, 1e-12)
        || !approx_abs(coe.argument_of_periapsis_rad, 0.0, 1e-12) || !approx_abs(coe.true_anomaly_rad, 0.0, 1e-12)) {
      spdlog::error("circular-equatorial reference mismatch");
      return 2;
    }
  }

  {
    const double a = 7200e3;
    const double e = 0.13;
    const double i = 0.9;
    const double raan = 1.2;
    const double argp = 0.7;
    const double nu = 2.1;
    const auto state = coe_to_state(a, e, i, raan, argp, nu, mu);
    const auto coe = state_eci_to_classical_orbital_elements(state.r_m, state.v_mps, mu);
    if (coe.status != Status::Ok) {
      spdlog::error("generic conversion failed");
      return 3;
    }
    if (!approx_abs(coe.semi_major_axis_m, a, 1e-3) || !approx_abs(coe.eccentricity, e, 1e-12)
        || !approx_angle(coe.inclination_rad, i, 1e-12) || !approx_angle(coe.raan_rad, raan, 1e-12)
        || !approx_angle(coe.argument_of_periapsis_rad, argp, 1e-12) || !approx_angle(coe.true_anomaly_rad, nu, 1e-12)) {
      spdlog::error("generic orbital-element mismatch");
      return 4;
    }

    const auto state_back = osculating_orbital_elements_to_state_eci(coe, mu);
    if (state_back.status != Status::Ok) {
      spdlog::error("generic reverse conversion failed");
      return 5;
    }
    if (!approx_abs(state_back.position_m.x, state.r_m.x, 1e-6) || !approx_abs(state_back.position_m.y, state.r_m.y, 1e-6)
        || !approx_abs(state_back.position_m.z, state.r_m.z, 1e-6) || !approx_abs(state_back.velocity_mps.x, state.v_mps.x, 1e-9)
        || !approx_abs(state_back.velocity_mps.y, state.v_mps.y, 1e-9) || !approx_abs(state_back.velocity_mps.z, state.v_mps.z, 1e-9)) {
      spdlog::error("generic state round-trip mismatch");
      return 6;
    }
  }

  {
    const StateVector bad_state{
        .position_m = Vec3{7000e3, 0.0, 0.0},
        .velocity_mps = Vec3{0.0, 7500.0, 0.0},
        .frame = Frame::ECEF,
    };
    const auto coe = state_to_classical_orbital_elements(bad_state, mu);
    if (coe.status != Status::InvalidInput) {
      spdlog::error("non-ECI frame should be rejected");
      return 7;
    }
  }

  {
    ClassicalOrbitalElements bad_elements{};
    bad_elements.semi_major_axis_m = 7000e3;
    bad_elements.eccentricity = 1.0;
    bad_elements.inclination_rad = 0.1;
    bad_elements.raan_rad = 0.2;
    bad_elements.argument_of_periapsis_rad = 0.3;
    bad_elements.true_anomaly_rad = 0.4;
    const auto state = osculating_orbital_elements_to_state_eci(bad_elements, mu);
    if (state.status != Status::InvalidInput) {
      spdlog::error("parabolic elements without p should be rejected");
      return 8;
    }
  }

  return 0;
}
