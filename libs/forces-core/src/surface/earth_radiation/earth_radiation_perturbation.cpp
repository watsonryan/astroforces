/**
 * @file earth_radiation_perturbation.cpp
 * @brief Earth Radiation perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/earth_radiation/earth_radiation_perturbation.hpp"

namespace astroforces::forces {

astroforces::forces::PerturbationContribution EarthRadiationPerturbationModel::evaluate(
    const astroforces::forces::PerturbationRequest& request) const {
  astroforces::forces::PerturbationContribution out{};
  out.name = name_;
  out.type = astroforces::forces::PerturbationType::EarthRadiation;

  const auto* sc = request.spacecraft ? request.spacecraft : default_spacecraft_;
  if (!sc) {
    out.status = astroforces::core::Status::InvalidInput;
    return out;
  }

  const auto earth_radiation = earth_radiation_.evaluate(request.state, *sc);
  out.acceleration_mps2 = earth_radiation.acceleration_mps2;
  out.status = earth_radiation.status;
  return out;
}

}  // namespace astroforces::forces
