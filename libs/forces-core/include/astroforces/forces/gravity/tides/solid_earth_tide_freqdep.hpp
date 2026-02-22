/**
 * @file solid_earth_tide_freqdep.hpp
 * @brief Frequency-dependent solid Earth tide (IERS2010 6.5 style) helpers.
 * @author Watosn
 */
#pragma once

#include <Eigen/Dense>

namespace astroforces::forces::tides {

void add_solid_earth_tide2_delta(double jd_utc, double gmst_rad, Eigen::MatrixXd& dC, Eigen::MatrixXd& dS);

}  // namespace astroforces::forces::tides
