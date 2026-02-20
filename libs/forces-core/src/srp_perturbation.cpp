/**
 * @file srp_perturbation.cpp
 * @brief SRP perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/srp/srp_perturbation.hpp"

namespace astroforces::srp {

astroforces::forces::PerturbationContribution SrpPerturbationModel::evaluate(
    const astroforces::forces::PerturbationRequest& request) const {
  astroforces::forces::PerturbationContribution out{};
  out.name = name_;
  out.type = astroforces::forces::PerturbationType::SRP;

  if (!srp_) {
    out.status = astroforces::atmo::Status::DataUnavailable;
    return out;
  }

  const auto* sc = request.spacecraft ? request.spacecraft : default_spacecraft_;
  if (!sc) {
    out.status = astroforces::atmo::Status::InvalidInput;
    return out;
  }

  const auto srp = srp_->evaluate(request.state, *sc);
  out.acceleration_mps2 = srp.acceleration_mps2;
  out.status = srp.status;
  return out;
}

}  // namespace astroforces::srp

