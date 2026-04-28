/**
 * @file solid_earth_tide.hpp
 * @brief Solid Earth tide contribution helpers for gravity SH models.
 * @author Watson
 */
#pragma once

#include <Eigen/Dense>

#include "astroforces/core/types.hpp"

namespace astroforces::forces::tides {

void add_solid_earth_tide1_delta(const astroforces::core::Vec3& body_ecef_m,
                                 double mu_body_m3_s2,
                                 double mu_earth_m3_s2,
                                 double radius_m,
                                 int max_deg,
                                 Eigen::MatrixXd& dC,
                                 Eigen::MatrixXd& dS);

}  // namespace astroforces::forces::tides
