/**
 * @file gravity_sph_perturbation.cpp
 * @brief Gravity perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/forces/gravity/gravity_sph_perturbation.hpp"

namespace astroforces::forces {

PerturbationContribution GravitySphPerturbationModel::evaluate(const PerturbationRequest& request) const {
  PerturbationContribution out{};
  out.name = name_;
  out.type = PerturbationType::Gravity;

  if (!gravity_) {
    out.status = astroforces::core::Status::DataUnavailable;
    return out;
  }

  const auto g = gravity_->evaluate(request.state);
  out.acceleration_mps2 = g.acceleration_mps2;
  out.status = g.status;
  return out;
}

}  // namespace astroforces::forces
