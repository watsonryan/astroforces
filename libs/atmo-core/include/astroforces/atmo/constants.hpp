/**
 * @file constants.hpp
 * @brief Shared physical constants for astrodynamics models.
 * @author Watosn
 */
#pragma once

namespace astroforces::atmo::constants {

inline constexpr double kEarthRadiusWgs84M = 6378137.0;
inline constexpr double kEarthEquatorialRadiusM = 6378136.3;
inline constexpr double kEarthMuM3S2 = 3.986004418e14;
inline constexpr double kSunMuM3S2 = 1.32712440018e20;
inline constexpr double kMoonMuM3S2 = 4.9048695e12;
inline constexpr double kEarthC20FullyNormalized = -4.84165143790815e-4;
inline constexpr double kAstronomicalUnitM = 149597870700.0;
inline constexpr double kSolarRadiationPressureAt1AuPa = 4.56e-6;
inline constexpr double kEarthRadiationPressureAtEarthRadiusPa = 8.0e-7;

}  // namespace astroforces::atmo::constants
