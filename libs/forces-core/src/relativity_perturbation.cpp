/**
 * @file relativity_perturbation.cpp
 * @brief Relativity perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/forces/relativity_perturbation.hpp"

namespace astroforces::forces {

PerturbationContribution RelativityPerturbationModel::evaluate(const PerturbationRequest& request) const {
  PerturbationContribution out{};
  out.name = name_;
  out.type = PerturbationType::Relativity;

  if (!relativity_) {
    out.status = astroforces::atmo::Status::DataUnavailable;
    return out;
  }

  const auto rel = relativity_->evaluate(request.state);
  out.acceleration_mps2 = rel.acceleration_mps2;
  out.status = rel.status;
  return out;
}

}  // namespace astroforces::forces

