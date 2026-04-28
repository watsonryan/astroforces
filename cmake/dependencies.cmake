# Author: Watson

include(${CMAKE_CURRENT_LIST_DIR}/CPM.cmake)

set(ASTROFORCES_NRLMSIS21_REPO "https://github.com/watsonryan/nrlmsis-2_1.git" CACHE STRING "NRLMSIS 2.1 repo")
set(ASTROFORCES_DTM2020_REPO "https://github.com/watsonryan/dtm2020.git" CACHE STRING "DTM2020 repo")
set(ASTROFORCES_HWM14_REPO "https://github.com/watsonryan/hwm14.git" CACHE STRING "HWM14 repo")
set(ASTROFORCES_JPL_EPH_REPO "https://github.com/watsonryan/jplEphem.git" CACHE STRING "jplEphem repo")
set(ASTROFORCES_EXTERNAL_MODELS_TAG "main" CACHE STRING "Branch/tag for external model repos")

set(ASTROFORCES_EIGEN_REPO "https://gitlab.com/libeigen/eigen.git" CACHE STRING "Eigen repository URL")
set(ASTROFORCES_EIGEN_TAG "5.0.0" CACHE STRING "Eigen git tag")
set(ASTROFORCES_CLI11_REPO "https://github.com/CLIUtils/CLI11.git" CACHE STRING "CLI11 repository URL")
set(ASTROFORCES_CLI11_TAG "v2.6.1" CACHE STRING "CLI11 git tag")

CPMAddPackage(
  NAME eigen
  GIT_REPOSITORY ${ASTROFORCES_EIGEN_REPO}
  GIT_TAG ${ASTROFORCES_EIGEN_TAG}
)
find_package(fmt CONFIG QUIET)
if(NOT TARGET fmt::fmt)
  CPMAddPackage(
    NAME fmt
    URL https://github.com/fmtlib/fmt/archive/refs/tags/12.1.0.tar.gz
    OPTIONS
      "FMT_DOC OFF"
      "FMT_TEST OFF"
      "FMT_FUZZ OFF"
  )
endif()
if(TARGET fmt AND NOT TARGET fmt::fmt)
  add_library(fmt::fmt ALIAS fmt)
endif()

CPMAddPackage(
  NAME cli11
  GIT_REPOSITORY ${ASTROFORCES_CLI11_REPO}
  GIT_TAG ${ASTROFORCES_CLI11_TAG}
  OPTIONS
    "CLI11_BUILD_TESTS OFF"
    "CLI11_BUILD_EXAMPLES OFF"
    "CLI11_BUILD_DOCS OFF"
)
CPMAddPackage(
  NAME spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.14.1
  OPTIONS
    "SPDLOG_FMT_EXTERNAL ON"
)
if(TARGET spdlog AND NOT TARGET spdlog::spdlog)
  add_library(spdlog::spdlog ALIAS spdlog)
endif()

CPMAddPackage(
  NAME nrlmsis21_ext
  GIT_REPOSITORY ${ASTROFORCES_NRLMSIS21_REPO}
  GIT_TAG ${ASTROFORCES_EXTERNAL_MODELS_TAG}
  OPTIONS
    "MSIS21_BUILD_TESTS OFF"
    "MSIS21_BUILD_CLI OFF"
)

CPMAddPackage(
  NAME dtm2020_ext
  GIT_REPOSITORY ${ASTROFORCES_DTM2020_REPO}
  GIT_TAG ${ASTROFORCES_EXTERNAL_MODELS_TAG}
  OPTIONS
    "DTM2020_BUILD_TESTS OFF"
    "DTM2020_ENABLE_RESEARCH OFF"
)

CPMAddPackage(
  NAME hwm14_ext
  GIT_REPOSITORY ${ASTROFORCES_HWM14_REPO}
  GIT_TAG ${ASTROFORCES_EXTERNAL_MODELS_TAG}
  OPTIONS
    "HWM14_BUILD_TESTS OFF"
    "HWM14_BUILD_EXAMPLES OFF"
)

CPMAddPackage(
  NAME jpl_eph_ext
  GIT_REPOSITORY ${ASTROFORCES_JPL_EPH_REPO}
  GIT_TAG ${ASTROFORCES_EXTERNAL_MODELS_TAG}
  OPTIONS
    "JPL_EPH_BUILD_TESTS OFF"
    "JPL_EPH_BUILD_TOOLS OFF"
    "JPL_EPH_FETCH_DEPS OFF"
)

set(ASTROFORCES_EIGEN_TARGET "")
if(TARGET Eigen3::Eigen)
  set(ASTROFORCES_EIGEN_TARGET Eigen3::Eigen)
elseif(TARGET eigen)
  set(ASTROFORCES_EIGEN_TARGET eigen)
endif()

set(ASTROFORCES_CLI11_TARGET "")
if(TARGET CLI11::CLI11)
  set(ASTROFORCES_CLI11_TARGET CLI11::CLI11)
endif()
