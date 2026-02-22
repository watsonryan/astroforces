/**
 * @file exponential_atmosphere.hpp
 * @brief Basic exponential atmosphere model for integration testing.
 * @author Watosn
 */
#pragma once

#include "astroforces/atmo/interfaces.hpp"

namespace astroforces::models {

class ExponentialAtmosphereModel final : public astroforces::core::IAtmosphereModel {
 public:
  ExponentialAtmosphereModel(double rho0_kg_m3, double h0_m, double scale_height_m, double temperature_k)
      : rho0_(rho0_kg_m3), h0_(h0_m), hs_(scale_height_m), t_(temperature_k) {}

  [[nodiscard]] astroforces::core::AtmosphereSample evaluate(const astroforces::core::StateVector& state,
                                                          const astroforces::core::WeatherIndices& weather) const override;

 private:
  double rho0_{};
  double h0_{};
  double hs_{};
  double t_{};
};

class ZeroWindModel final : public astroforces::core::IWindModel {
 public:
  [[nodiscard]] astroforces::core::WindSample evaluate(const astroforces::core::StateVector& state,
                                                    const astroforces::core::WeatherIndices& weather) const override;
};

}  // namespace astroforces::models
