/**
 * @file perturbation.cpp
 * @brief Generic perturbation model aggregation implementation.
 * @author Watosn
 */

#include "astroforces/forces/perturbation.hpp"

namespace astroforces::forces {

void PerturbationStack::add(std::unique_ptr<IPerturbationModel> model) {
  if (!model) {
    return;
  }
  models_.push_back(std::move(model));
}

PerturbationResult PerturbationStack::evaluate(const PerturbationRequest& request) const {
  PerturbationResult out{};
  for (const auto& model : models_) {
    const auto c = model->evaluate(request);
    out.contributions.push_back(c);
    if (c.status != astroforces::atmo::Status::Ok && out.status == astroforces::atmo::Status::Ok) {
      out.status = c.status;
    }
    out.total_acceleration_mps2 = out.total_acceleration_mps2 + c.acceleration_mps2;
  }
  return out;
}

}  // namespace astroforces::forces
