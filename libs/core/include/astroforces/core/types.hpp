/**
 * @file types.hpp
 * @brief Core domain types for astroforces.
 * @author Watosn
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

namespace astroforces::core {

/**
 * @brief Supported coordinate frames for states/vectors.
 */
enum class Frame : std::uint8_t { ECI, ECEF, NED, BODY };

/**
 * @brief Standard status code used by model outputs.
 */
enum class Status : std::uint8_t { Ok, InvalidInput, NotImplemented, DataUnavailable, NumericalError };

/**
 * @brief Source tag for weather data provenance.
 */
enum class WeatherSource : std::uint8_t { Unknown, StaticProvider, CelesTrakLast5YearsCsv };

/**
 * @brief Cartesian 3-vector.
 */
struct Vec3 {
  double x{};
  double y{};
  double z{};
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3{a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3{a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, const Vec3& v) { return Vec3{s * v.x, s * v.y, s * v.z}; }
inline Vec3 operator*(const Vec3& v, double s) { return s * v; }
inline Vec3 operator/(const Vec3& v, double s) { return Vec3{v.x / s, v.y / s, v.z / s}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3& v) { return std::sqrt(dot(v, v)); }

/**
 * @brief UTC epoch expressed as seconds since Unix epoch.
 */
struct Epoch {
  double utc_seconds{};
};

/**
 * @brief Full spacecraft state at an epoch.
 *
 * `body_from_frame_dcm` is row-major and transforms vectors from `frame` into BODY.
 */
struct StateVector {
  Epoch epoch{};
  Vec3 position_m{};
  Vec3 velocity_mps{};
  Frame frame{Frame::ECI};
  std::array<double, 9> body_from_frame_dcm{1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
};

/**
 * @brief Geodetic-like point used for atmospheric/wind model inputs.
 */
struct GeodeticPoint {
  double lat_deg{};
  double lon_deg{};
  double alt_m{};
};

/**
 * @brief Earth orientation parameters at an epoch.
 */
struct EarthOrientation {
  double xp_rad{};
  double yp_rad{};
  double dut1_s{};
  double lod_s{};
  double dX_rad{};
  double dY_rad{};
};

/**
 * @brief Time derivative of Earth orientation parameters.
 */
struct EarthOrientationRate {
  double xp_rad_s{};
  double yp_rad_s{};
  double dut1_s_s{};
  double lod_s_s{};
  double dX_rad_s{};
  double dY_rad_s{};
};

/**
 * @brief Celestial Intermediate Pole state (X,Y,s).
 */
struct CelestialIntermediatePole {
  double x_rad{};
  double y_rad{};
  double s_rad{};
};

/**
 * @brief Time derivative of CIP state (Xdot,Ydot,sdot).
 */
struct CelestialIntermediatePoleRate {
  double x_rad_s{};
  double y_rad_s{};
  double s_rad_s{};
};

/**
 * @brief Space-weather forcing bundle shared across atmospheric models.
 */
struct WeatherIndices {
  double f107{};
  double f107a{};
  double ap{};
  double kp{};
  double ap_3h_current{};
  double kp_3h_current{};
  std::array<double, 8> ap_3h_utc{};
  std::array<double, 8> kp_3h_utc{};
  std::array<double, 7> ap_msis_history{};
  bool has_ap_msis_history{};
  bool f107_observed{true};
  bool geomagnetic_observed{true};
  WeatherSource source{WeatherSource::Unknown};
  bool interpolated{};
  bool extrapolated{};
  Status status{Status::Ok};
};

/**
 * @brief Atmospheric density/temperature sample.
 */
struct AtmosphereSample {
  double density_kg_m3{};
  double temperature_k{};
  Status status{Status::Ok};
};

/**
 * @brief Neutral wind sample.
 */
struct WindSample {
  Vec3 velocity_mps{};
  Frame frame{Frame::ECI};
  Status status{Status::Ok};
};

}  // namespace astroforces::core
