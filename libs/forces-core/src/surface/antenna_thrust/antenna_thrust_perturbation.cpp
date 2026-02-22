/**
 * @file antenna_thrust_perturbation.cpp
 * @brief Antenna thrust perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_perturbation.hpp"

namespace astroforces::forces {

astroforces::forces::PerturbationContribution AntennaThrustPerturbationModel::evaluate(
    const astroforces::forces::PerturbationRequest& request) const {
  astroforces::forces::PerturbationContribution out{};
  out.name = name_;
  out.type = astroforces::forces::PerturbationType::AntennaThrust;

  const auto* sc = request.spacecraft ? request.spacecraft : default_spacecraft_;
  if (!sc) {
    out.status = astroforces::core::Status::InvalidInput;
    return out;
  }

  const auto antenna = antenna_thrust_.evaluate(request.state, *sc);
  out.acceleration_mps2 = antenna.acceleration_mps2;
  out.status = antenna.status;
  return out;
}

}  // namespace astroforces::forces

