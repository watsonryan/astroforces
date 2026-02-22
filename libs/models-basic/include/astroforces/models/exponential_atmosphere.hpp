/**
 * @file exponential_atmosphere.hpp
 * @brief Basic exponential atmosphere model for integration testing.
 * @author Watosn
 */
#pragma once

#include "astroforces/core/interfaces.hpp"

namespace astroforces::models {

/**
 * @brief Simple exponential atmosphere model for testing/smoke workflows.
 */
class ExponentialAtmosphereModel final : public astroforces::core::IAtmosphereModel {
 public:
  /**
   * @brief Construct model with reference density/height and scale height.
   */
  ExponentialAtmosphereModel(double rho0_kg_m3, double h0_m, double scale_height_m, double temperature_k)
      : rho0_(rho0_kg_m3), h0_(h0_m), hs_(scale_height_m), t_(temperature_k) {}

  /**
   * @brief Evaluate exponential density and fixed temperature.
   */
  [[nodiscard]] astroforces::core::AtmosphereSample evaluate(const astroforces::core::StateVector& state,
                                                          const astroforces::core::WeatherIndices& weather) const override;

 private:
  double rho0_{};
  double h0_{};
  double hs_{};
  double t_{};
};

/**
 * @brief Zero-wind model returning a null velocity vector.
 */
class ZeroWindModel final : public astroforces::core::IWindModel {
 public:
  /**
   * @brief Evaluate zero wind.
   */
  [[nodiscard]] astroforces::core::WindSample evaluate(const astroforces::core::StateVector& state,
                                                    const astroforces::core::WeatherIndices& weather) const override;
};

}  // namespace astroforces::models
