/**
 * @file orbital_elements.hpp
 * @brief Classical orbital-element utilities.
 * @author Watosn
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/math_utils.hpp"
#include "astroforces/core/types.hpp"

namespace astroforces::core {

/**
 * @brief Classical osculating orbital elements derived from an instantaneous Cartesian state.
 *
 * The returned element set is `(a, e, i, RAAN, omega, nu)`, with all angles in radians.
 *
 * Deterministic conventions are applied at the usual singularities:
 * - equatorial orbit: `raan_rad = 0`
 * - circular inclined orbit: `argument_of_periapsis_rad = 0`, `true_anomaly_rad = argument of latitude`
 * - eccentric equatorial orbit: `raan_rad = 0`, `argument_of_periapsis_rad = longitude of periapsis`
 * - circular equatorial orbit: `raan_rad = 0`, `argument_of_periapsis_rad = 0`,
 *   `true_anomaly_rad = true longitude`
 */
struct ClassicalOrbitalElements {
  double semi_major_axis_m{};
  double semi_latus_rectum_m{};
  double eccentricity{};
  double inclination_rad{};
  double raan_rad{};
  double argument_of_periapsis_rad{};
  double true_anomaly_rad{};
  double eccentric_anomaly_rad{};
  double mean_anomaly_rad{};
  bool circular{};
  bool equatorial{};
  Status status{Status::Ok};
};

/**
 * @brief ECI Cartesian state returned by orbital-element conversion utilities.
 */
struct CartesianOrbitState {
  Vec3 position_m{};
  Vec3 velocity_mps{};
  Frame frame{Frame::ECI};
  Status status{Status::Ok};
};

namespace detail {

inline double clamp_unit(const double value) {
  return std::clamp(value, -1.0, 1.0);
}

inline double wrap_to_2pi(double angle_rad) {
  angle_rad = std::fmod(angle_rad, constants::kTwoPi);
  if (angle_rad < 0.0) {
    angle_rad += constants::kTwoPi;
  }
  return angle_rad;
}

}  // namespace detail

/**
 * @brief Convert an ECI Cartesian state to classical osculating orbital elements.
 *
 * This function assumes the input state is inertial and centered on the body associated with `mu_m3_s2`.
 * For Earth-centered use, the default gravitational parameter is used.
 *
 * @param position_m ECI position vector in meters.
 * @param velocity_mps ECI velocity vector in meters per second.
 * @param mu_m3_s2 Central-body gravitational parameter in m^3/s^2.
 * @return Classical orbital elements with `status` set.
 */
[[nodiscard]] inline ClassicalOrbitalElements state_eci_to_classical_orbital_elements(const Vec3& position_m,
                                                                                       const Vec3& velocity_mps,
                                                                                       const double mu_m3_s2 = constants::kEarthMuM3S2) {
  constexpr double kParabolicTol = 1e-12;

  ClassicalOrbitalElements out{};
  if (!(mu_m3_s2 > 0.0)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const double r_mag = norm(position_m);
  const double v_mag = norm(velocity_mps);
  if (!(r_mag > 0.0) || !std::isfinite(r_mag) || !std::isfinite(v_mag)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const Vec3 h = vec_cross(position_m, velocity_mps);
  const double h_mag = norm(h);
  if (!(h_mag > 0.0) || !std::isfinite(h_mag)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const Vec3 k_hat{0.0, 0.0, 1.0};
  const Vec3 node = vec_cross(k_hat, h);
  const double node_mag = norm(node);

  const Vec3 e_vec = vec_cross(velocity_mps, h) / mu_m3_s2 - position_m / r_mag;
  const double e = norm(e_vec);
  const double specific_energy = 0.5 * v_mag * v_mag - mu_m3_s2 / r_mag;
  const double p = h_mag * h_mag / mu_m3_s2;

  out.semi_latus_rectum_m = p;
  out.eccentricity = e;
  out.circular = e < 1e-10;
  out.equatorial = node_mag < 1e-10;

  const double cos_i = detail::clamp_unit(h.z / h_mag);
  out.inclination_rad = std::acos(cos_i);

  if (std::abs(specific_energy) < kParabolicTol) {
    out.semi_major_axis_m = std::numeric_limits<double>::infinity();
  } else {
    out.semi_major_axis_m = -mu_m3_s2 / (2.0 * specific_energy);
  }

  if (out.equatorial) {
    out.raan_rad = 0.0;
  } else {
    out.raan_rad = detail::wrap_to_2pi(std::atan2(node.y, node.x));
  }

  if (!out.circular && !out.equatorial) {
    const double cos_argp = detail::clamp_unit(dot(node, e_vec) / (node_mag * e));
    double argp = std::acos(cos_argp);
    if (e_vec.z < 0.0) {
      argp = constants::kTwoPi - argp;
    }
    out.argument_of_periapsis_rad = detail::wrap_to_2pi(argp);
  } else if (!out.circular && out.equatorial) {
    out.argument_of_periapsis_rad = detail::wrap_to_2pi(std::atan2(e_vec.y, e_vec.x));
  } else {
    out.argument_of_periapsis_rad = 0.0;
  }

  if (!out.circular) {
    const double cos_nu = detail::clamp_unit(dot(e_vec, position_m) / (e * r_mag));
    double nu = std::acos(cos_nu);
    if (dot(position_m, velocity_mps) < 0.0) {
      nu = constants::kTwoPi - nu;
    }
    out.true_anomaly_rad = detail::wrap_to_2pi(nu);
  } else if (!out.equatorial) {
    const double cos_u = detail::clamp_unit(dot(node, position_m) / (node_mag * r_mag));
    double u = std::acos(cos_u);
    if (position_m.z < 0.0) {
      u = constants::kTwoPi - u;
    }
    out.true_anomaly_rad = detail::wrap_to_2pi(u);
  } else {
    out.true_anomaly_rad = detail::wrap_to_2pi(std::atan2(position_m.y, position_m.x));
  }

  if (e < 1.0 - 1e-12) {
    const double sin_half_nu = std::sin(out.true_anomaly_rad / 2.0);
    const double cos_half_nu = std::cos(out.true_anomaly_rad / 2.0);
    out.eccentric_anomaly_rad = 2.0 * std::atan2(std::sqrt(1.0 - e) * sin_half_nu,
                                                 std::sqrt(1.0 + e) * cos_half_nu);
    out.eccentric_anomaly_rad = detail::wrap_to_2pi(out.eccentric_anomaly_rad);
    out.mean_anomaly_rad =
        detail::wrap_to_2pi(out.eccentric_anomaly_rad - e * std::sin(out.eccentric_anomaly_rad));
  } else {
    out.eccentric_anomaly_rad = std::numeric_limits<double>::quiet_NaN();
    out.mean_anomaly_rad = std::numeric_limits<double>::quiet_NaN();
  }

  if (!std::isfinite(out.semi_latus_rectum_m) || !std::isfinite(out.eccentricity)
      || !std::isfinite(out.inclination_rad) || !std::isfinite(out.raan_rad)
      || !std::isfinite(out.argument_of_periapsis_rad) || !std::isfinite(out.true_anomaly_rad)) {
    out.status = Status::NumericalError;
    return out;
  }

  out.status = Status::Ok;
  return out;
}

/**
 * @brief Alias for converting an ECI state to standard classical osculating elements.
 */
[[nodiscard]] inline ClassicalOrbitalElements state_eci_to_osculating_orbital_elements(const Vec3& position_m,
                                                                                        const Vec3& velocity_mps,
                                                                                        const double mu_m3_s2 = constants::kEarthMuM3S2) {
  return state_eci_to_classical_orbital_elements(position_m, velocity_mps, mu_m3_s2);
}

/**
 * @brief Convert classical osculating orbital elements to an ECI Cartesian state.
 *
 * Uses the standard perifocal-to-inertial rotation based on `(RAAN, omega, i)`.
 * The input is interpreted as osculating elements about the central body associated with `mu_m3_s2`.
 *
 * `semi_latus_rectum_m` is preferred when provided. Otherwise, for non-parabolic cases it is
 * reconstructed from `a * (1 - e^2)`.
 *
 * @param elements Classical osculating orbital elements.
 * @param mu_m3_s2 Central-body gravitational parameter in m^3/s^2.
 * @return ECI Cartesian state with `status` set.
 */
[[nodiscard]] inline CartesianOrbitState classical_orbital_elements_to_state_eci(const ClassicalOrbitalElements& elements,
                                                                                  const double mu_m3_s2 = constants::kEarthMuM3S2) {
  CartesianOrbitState out{};
  if (!(mu_m3_s2 > 0.0) || !std::isfinite(elements.eccentricity) || !std::isfinite(elements.inclination_rad)
      || !std::isfinite(elements.raan_rad) || !std::isfinite(elements.argument_of_periapsis_rad)
      || !std::isfinite(elements.true_anomaly_rad)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const double e = elements.eccentricity;
  if (e < 0.0) {
    out.status = Status::InvalidInput;
    return out;
  }

  double p = elements.semi_latus_rectum_m;
  if (!(p > 0.0) || !std::isfinite(p)) {
    const double one_minus_e2 = 1.0 - e * e;
    if (std::abs(one_minus_e2) < 1e-12 || !std::isfinite(elements.semi_major_axis_m)) {
      out.status = Status::InvalidInput;
      return out;
    }
    p = elements.semi_major_axis_m * one_minus_e2;
  }

  if (!(p > 0.0) || !std::isfinite(p)) {
    out.status = Status::InvalidInput;
    return out;
  }

  const double nu = elements.true_anomaly_rad;
  const double cos_nu = std::cos(nu);
  const double sin_nu = std::sin(nu);
  const double denom = 1.0 + e * cos_nu;
  if (std::abs(denom) < 1e-14 || !std::isfinite(denom)) {
    out.status = Status::NumericalError;
    return out;
  }

  const double r = p / denom;
  const double h = std::sqrt(mu_m3_s2 * p);
  if (!(r > 0.0) || !std::isfinite(r) || !(h > 0.0) || !std::isfinite(h)) {
    out.status = Status::NumericalError;
    return out;
  }

  const Vec3 r_pf{r * cos_nu, r * sin_nu, 0.0};
  const Vec3 v_pf{-mu_m3_s2 / h * sin_nu, mu_m3_s2 / h * (e + cos_nu), 0.0};

  const double cO = std::cos(elements.raan_rad);
  const double sO = std::sin(elements.raan_rad);
  const double ci = std::cos(elements.inclination_rad);
  const double si = std::sin(elements.inclination_rad);
  const double cw = std::cos(elements.argument_of_periapsis_rad);
  const double sw = std::sin(elements.argument_of_periapsis_rad);

  const double q11 = cO * cw - sO * sw * ci;
  const double q12 = -cO * sw - sO * cw * ci;
  const double q21 = sO * cw + cO * sw * ci;
  const double q22 = -sO * sw + cO * cw * ci;
  const double q31 = sw * si;
  const double q32 = cw * si;

  out.position_m = Vec3{
      q11 * r_pf.x + q12 * r_pf.y,
      q21 * r_pf.x + q22 * r_pf.y,
      q31 * r_pf.x + q32 * r_pf.y,
  };
  out.velocity_mps = Vec3{
      q11 * v_pf.x + q12 * v_pf.y,
      q21 * v_pf.x + q22 * v_pf.y,
      q31 * v_pf.x + q32 * v_pf.y,
  };

  if (!std::isfinite(out.position_m.x) || !std::isfinite(out.position_m.y) || !std::isfinite(out.position_m.z)
      || !std::isfinite(out.velocity_mps.x) || !std::isfinite(out.velocity_mps.y)
      || !std::isfinite(out.velocity_mps.z)) {
    out.status = Status::NumericalError;
    return out;
  }

  out.status = Status::Ok;
  return out;
}

/**
 * @brief Alias for converting standard osculating elements to an ECI Cartesian state.
 */
[[nodiscard]] inline CartesianOrbitState osculating_orbital_elements_to_state_eci(const ClassicalOrbitalElements& elements,
                                                                                   const double mu_m3_s2 = constants::kEarthMuM3S2) {
  return classical_orbital_elements_to_state_eci(elements, mu_m3_s2);
}

/**
 * @brief Convert a full state vector to classical osculating orbital elements.
 *
 * Requires `state.frame == Frame::ECI`; otherwise `Status::InvalidInput` is returned.
 *
 * @param state Input state vector in ECI.
 * @param mu_m3_s2 Central-body gravitational parameter in m^3/s^2.
 * @return Classical orbital elements with `status` set.
 */
[[nodiscard]] inline ClassicalOrbitalElements state_to_classical_orbital_elements(const StateVector& state,
                                                                                   const double mu_m3_s2 = constants::kEarthMuM3S2) {
  if (state.frame != Frame::ECI) {
    return ClassicalOrbitalElements{.status = Status::InvalidInput};
  }
  return state_eci_to_classical_orbital_elements(state.position_m, state.velocity_mps, mu_m3_s2);
}

/**
 * @brief Alias for converting a state vector to standard classical osculating elements.
 *
 * Requires `state.frame == Frame::ECI`; otherwise `Status::InvalidInput` is returned.
 */
[[nodiscard]] inline ClassicalOrbitalElements state_to_osculating_orbital_elements(const StateVector& state,
                                                                                    const double mu_m3_s2 = constants::kEarthMuM3S2) {
  return state_to_classical_orbital_elements(state, mu_m3_s2);
}

}  // namespace astroforces::core
