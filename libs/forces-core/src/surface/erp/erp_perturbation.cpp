/**
 * @file erp_perturbation.cpp
 * @brief ERP perturbation wrapper implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/erp/erp_perturbation.hpp"

namespace astroforces::forces {

astroforces::forces::PerturbationContribution ErpPerturbationModel::evaluate(
    const astroforces::forces::PerturbationRequest& request) const {
  astroforces::forces::PerturbationContribution out{};
  out.name = name_;
  out.type = astroforces::forces::PerturbationType::ERP;

  const auto* sc = request.spacecraft ? request.spacecraft : default_spacecraft_;
  if (!sc) {
    out.status = astroforces::core::Status::InvalidInput;
    return out;
  }

  const auto erp = erp_.evaluate(request.state, *sc);
  out.acceleration_mps2 = erp.acceleration_mps2;
  out.status = erp.status;
  return out;
}

}  // namespace astroforces::forces

