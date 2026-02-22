/**
 * @file zero_wind.cpp
 * @brief Zero neutral wind model implementation.
 * @author Watosn
 */

#include "astroforces/models/exponential_atmosphere.hpp"

namespace astroforces::models {

astroforces::core::WindSample ZeroWindModel::evaluate(const astroforces::core::StateVector& state,
                                                   const astroforces::core::WeatherIndices& /*weather*/) const {
  return astroforces::core::WindSample{
      .velocity_mps = astroforces::core::Vec3{},
      .frame = state.frame,
      .status = astroforces::core::Status::Ok};
}

}  // namespace astroforces::models
