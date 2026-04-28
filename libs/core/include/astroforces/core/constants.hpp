/**
 * @file constants.hpp
 * @brief Shared physical constants for astrodynamics models.
 * @author Watson
 */
#pragma once

namespace astroforces::core::constants {

inline constexpr double kPi = 3.1415926535897932384626433832795;
inline constexpr double kTwoPi = 2.0 * kPi;
inline constexpr double kDegToRad = kPi / 180.0;
inline constexpr double kArcsecToRad = kDegToRad / 3600.0;
inline constexpr double kSecondsPerDay = 86400.0;
inline constexpr double kJ2000Jd = 2451545.0;

inline constexpr double kEarthRadiusWgs84M = 6378137.0;
inline constexpr double kEarthEquatorialRadiusM = 6378136.3;
inline constexpr double kSunRadiusM = 6.9634e8;
inline constexpr double kMoonRadiusM = 1.7374e6;
inline constexpr double kEarthMuM3S2 = 3.986004418e14;
inline constexpr double kEarthRotationRateRadPerSec = 7.2921150e-5;
inline constexpr double kEarthJ2 = 1.08262668e-3;
inline constexpr double kSunMuM3S2 = 1.32712440018e20;
inline constexpr double kMoonMuM3S2 = 4.9048695e12;
inline constexpr double kEarthC20FullyNormalized = -4.84165143790815e-4;
inline constexpr double kAstronomicalUnitM = 149597870700.0;
inline constexpr double kSpeedOfLightMps = 299792458.0;
inline constexpr double kEarthLoveNumberK2 = 0.3;
inline constexpr double kSolarRadiationPressureAt1AuPa = 4.56e-6;
inline constexpr double kSolarIrradianceAt1AuWm2 = 1361.0;
inline constexpr double kEarthBondAlbedo = 0.3;
inline constexpr double kEarthIrFluxWm2 = 237.0;
inline constexpr double kEarthRadiationPressureAtEarthRadiusPa = 8.0e-7;

}  // namespace astroforces::core::constants
