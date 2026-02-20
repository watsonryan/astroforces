/**
 * @file drag_perturbation.cpp
 * @brief Drag perturbation model implementation.
 * @author Watosn
 */

#include "astroforces/drag/drag_perturbation.hpp"

namespace astroforces::drag {

astroforces::forces::PerturbationContribution DragPerturbationModel::evaluate(
    const astroforces::forces::PerturbationRequest& request) const {
  astroforces::forces::PerturbationContribution out{};
  out.name = name_;
  out.type = astroforces::forces::PerturbationType::Drag;

  const auto* sc = request.spacecraft ? request.spacecraft : default_spacecraft_;
  if (!sc) {
    out.status = astroforces::atmo::Status::InvalidInput;
    return out;
  }

  const auto drag = drag_.evaluate(request.state, *sc);
  out.acceleration_mps2 = drag.acceleration_mps2;
  out.status = drag.status;
  return out;
}

}  // namespace astroforces::drag
