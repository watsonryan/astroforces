/**
 * @file exponential_atmosphere.cpp
 * @brief Basic atmosphere model implementation.
 * @author Watosn
 */

#include "astroforces/models/exponential_atmosphere.hpp"

#include <cmath>

#include "astroforces/atmo/constants.hpp"

namespace astroforces::models {

astroforces::core::AtmosphereSample ExponentialAtmosphereModel::evaluate(const astroforces::core::StateVector& state,
                                                                      const astroforces::core::WeatherIndices& /*weather*/) const {
  const double r = astroforces::core::norm(state.position_m);
  const double alt = r - astroforces::core::constants::kEarthRadiusWgs84M;
  const double rho = rho0_ * std::exp(-(alt - h0_) / hs_);
  return astroforces::core::AtmosphereSample{.density_kg_m3 = rho, .temperature_k = t_, .status = astroforces::core::Status::Ok};
}

}  // namespace astroforces::models
