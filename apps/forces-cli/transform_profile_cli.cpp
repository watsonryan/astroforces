/**
 * @file transform_profile_cli.cpp
 * @brief Micro-benchmark for rigorous GCRF<->ITRF transforms.
 * @author Watosn
 */

#include <chrono>
#include <cstdint>
#include <cstdlib>

#include <fmt/format.h>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/transforms.hpp"

namespace {

double parse_or(const char* s, const double fallback) {
  if (s == nullptr) {
    return fallback;
  }
  return std::strtod(s, nullptr);
}

std::uint64_t parse_or_u64(const char* s, const std::uint64_t fallback) {
  if (s == nullptr) {
    return fallback;
  }
  return static_cast<std::uint64_t>(std::strtoull(s, nullptr, 10));
}

}  // namespace

int main(int argc, char** argv) {
  const std::uint64_t iters = (argc >= 2) ? parse_or_u64(argv[1], 2000000ULL) : 2000000ULL;
  const double utc_seconds = (argc >= 3) ? parse_or(argv[2], 1000000000.0) : 1000000000.0;

  const double jd_utc = astroforces::core::utc_seconds_to_julian_date_utc(utc_seconds);
  const double jd_tt = astroforces::core::utc_seconds_to_julian_date_tt(utc_seconds);

  const astroforces::core::CelestialIntermediatePole cip{
      .x_rad = 0.2 * astroforces::core::constants::kArcsecToRad,
      .y_rad = -0.15 * astroforces::core::constants::kArcsecToRad,
      .s_rad = -0.01 * astroforces::core::constants::kArcsecToRad,
  };
  const astroforces::core::CelestialIntermediatePoleRate cip_rate{
      .x_rad_s = 1.0e-16,
      .y_rad_s = -1.0e-16,
      .s_rad_s = 5.0e-17,
  };
  const astroforces::core::EarthOrientation eop{
      .xp_rad = 0.12 * astroforces::core::constants::kArcsecToRad,
      .yp_rad = 0.19 * astroforces::core::constants::kArcsecToRad,
      .dut1_s = 0.083,
      .lod_s = 0.0009,
      .dX_rad = 0.01 * astroforces::core::constants::kArcsecToRad,
      .dY_rad = -0.01 * astroforces::core::constants::kArcsecToRad,
  };
  const astroforces::core::EarthOrientationRate eop_rate{
      .xp_rad_s = 2.0e-16,
      .yp_rad_s = -2.0e-16,
      .dut1_s_s = 0.0,
      .lod_s_s = 0.0,
      .dX_rad_s = 1.0e-16,
      .dY_rad_s = -1.0e-16,
  };

  const astroforces::core::Vec3 r_gcrf{7000e3, -1200e3, 2500e3};
  const astroforces::core::Vec3 v_gcrf{1200.0, 6800.0, -900.0};

  volatile double checksum_rotation = 0.0;
  volatile double checksum_position = 0.0;
  volatile double checksum_velocity = 0.0;
  volatile double checksum_position_cached = 0.0;
  volatile double checksum_velocity_cached = 0.0;

  {
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) {
      const auto rd =
          astroforces::core::gcrf_to_itrf_rotation_with_derivative_exact(jd_utc, jd_tt, cip, cip_rate, eop, eop_rate);
      checksum_rotation += rd.r(0, 0) + rd.r(1, 1) + rd.r(2, 2) + rd.dr(0, 0) + rd.dr(1, 1) + rd.dr(2, 2);
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    fmt::print("bench=rotation_with_derivative_exact iters={} total_ns={:.0f} ns_per_call={:.3f}\n",
               iters,
               ns,
               ns / static_cast<double>(iters));
  }

  {
    const auto ctx = astroforces::core::gcrf_to_itrf_transform_context(jd_utc, jd_tt, cip, eop);
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) {
      const astroforces::core::Vec3 r{
          r_gcrf.x + static_cast<double>(i & 7ULL),
          r_gcrf.y - static_cast<double>((i >> 3U) & 7ULL),
          r_gcrf.z + static_cast<double>((i >> 6U) & 7ULL),
      };
      const auto r_itrf = astroforces::core::gcrf_to_itrf_position(r, ctx);
      const auto r_back = astroforces::core::itrf_to_gcrf_position(r_itrf, ctx);
      checksum_position_cached += r_itrf.x + r_itrf.y + r_itrf.z + r_back.x + r_back.y + r_back.z;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    fmt::print("bench=position_roundtrip_cached iters={} total_ns={:.0f} ns_per_call={:.3f}\n",
               iters,
               ns,
               ns / static_cast<double>(iters));
  }

  {
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) {
      const astroforces::core::Vec3 r{
          r_gcrf.x + static_cast<double>(i & 7ULL),
          r_gcrf.y - static_cast<double>((i >> 3U) & 7ULL),
          r_gcrf.z + static_cast<double>((i >> 6U) & 7ULL),
      };
      const auto r_itrf = astroforces::core::gcrf_to_itrf_position(r, jd_utc, jd_tt, cip, eop);
      const auto r_back = astroforces::core::itrf_to_gcrf_position(r_itrf, jd_utc, jd_tt, cip, eop);
      checksum_position += r_itrf.x + r_itrf.y + r_itrf.z + r_back.x + r_back.y + r_back.z;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    fmt::print("bench=position_roundtrip iters={} total_ns={:.0f} ns_per_call={:.3f}\n",
               iters,
               ns,
               ns / static_cast<double>(iters));
  }

  {
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) {
      const astroforces::core::Vec3 r{
          r_gcrf.x + static_cast<double>(i & 7ULL),
          r_gcrf.y - static_cast<double>((i >> 3U) & 7ULL),
          r_gcrf.z + static_cast<double>((i >> 6U) & 7ULL),
      };
      const astroforces::core::Vec3 v{
          v_gcrf.x + static_cast<double>(i & 3ULL),
          v_gcrf.y - static_cast<double>((i >> 2U) & 3ULL),
          v_gcrf.z + static_cast<double>((i >> 4U) & 3ULL),
      };
      const auto v_itrf = astroforces::core::gcrf_to_itrf_velocity(r, v, jd_utc, jd_tt, cip, eop);
      const auto v_back = astroforces::core::itrf_to_gcrf_velocity(r, v_itrf, jd_utc, jd_tt, cip, eop);
      checksum_velocity += v_itrf.x + v_itrf.y + v_itrf.z + v_back.x + v_back.y + v_back.z;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    fmt::print("bench=velocity_roundtrip iters={} total_ns={:.0f} ns_per_call={:.3f}\n",
               iters,
               ns,
               ns / static_cast<double>(iters));
  }

  {
    const auto ctx = astroforces::core::gcrf_to_itrf_transform_context(jd_utc, jd_tt, cip, eop);
    const auto t0 = std::chrono::steady_clock::now();
    for (std::uint64_t i = 0; i < iters; ++i) {
      const astroforces::core::Vec3 r{
          r_gcrf.x + static_cast<double>(i & 7ULL),
          r_gcrf.y - static_cast<double>((i >> 3U) & 7ULL),
          r_gcrf.z + static_cast<double>((i >> 6U) & 7ULL),
      };
      const astroforces::core::Vec3 v{
          v_gcrf.x + static_cast<double>(i & 3ULL),
          v_gcrf.y - static_cast<double>((i >> 2U) & 3ULL),
          v_gcrf.z + static_cast<double>((i >> 4U) & 3ULL),
      };
      const auto v_itrf = astroforces::core::gcrf_to_itrf_velocity(r, v, ctx);
      const auto v_back = astroforces::core::itrf_to_gcrf_velocity(r, v_itrf, ctx);
      checksum_velocity_cached += v_itrf.x + v_itrf.y + v_itrf.z + v_back.x + v_back.y + v_back.z;
    }
    const auto t1 = std::chrono::steady_clock::now();
    const double ns = static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count());
    fmt::print("bench=velocity_roundtrip_cached iters={} total_ns={:.0f} ns_per_call={:.3f}\n",
               iters,
               ns,
               ns / static_cast<double>(iters));
  }

  fmt::print("checksum_rotation={:.17e}\n", static_cast<double>(checksum_rotation));
  fmt::print("checksum_position={:.17e}\n", static_cast<double>(checksum_position));
  fmt::print("checksum_velocity={:.17e}\n", static_cast<double>(checksum_velocity));
  fmt::print("checksum_position_cached={:.17e}\n", static_cast<double>(checksum_position_cached));
  fmt::print("checksum_velocity_cached={:.17e}\n", static_cast<double>(checksum_velocity_cached));
  return 0;
}
