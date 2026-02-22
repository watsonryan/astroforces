# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(astroforces VERSION 0.1.0 LANGUAGES CXX)

option(ASTROFORCES_BUILD_TESTS "Build tests" ON)
option(ASTROFORCES_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)
option(ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA "Fail tests when required external ephemeris data is missing" OFF)

set(ASTROFORCES_NRLMSIS21_REPO "https://github.com/watsonryan/nrlmsis-2_1.git" CACHE STRING "NRLMSIS 2.1 repo")
set(ASTROFORCES_DTM2020_REPO "https://github.com/watsonryan/dtm2020.git" CACHE STRING "DTM2020 repo")
set(ASTROFORCES_HWM14_REPO "https://github.com/watsonryan/hwm14.git" CACHE STRING "HWM14 repo")
set(ASTROFORCES_JPL_EPH_REPO "https://github.com/watsonryan/jplEphem.git" CACHE STRING "jplEphem repo")
set(ASTROFORCES_EXTERNAL_MODELS_TAG "main" CACHE STRING "Branch/tag for external model repos")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_SOURCE_DIR}/cmake/CPM.cmake)
set(ASTROFORCES_EIGEN_REPO "https://gitlab.com/libeigen/eigen.git" CACHE STRING "Eigen repository URL")
set(ASTROFORCES_EIGEN_TAG "5.0.0" CACHE STRING "Eigen git tag")

CPMAddPackage(
  NAME eigen
  GIT_REPOSITORY ${ASTROFORCES_EIGEN_REPO}
  GIT_TAG ${ASTROFORCES_EIGEN_TAG}
)
CPMAddPackage(
  NAME fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 10.2.1
)
CPMAddPackage(
  NAME spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.14.1
  OPTIONS
    "SPDLOG_FMT_EXTERNAL ON"
)

set(ASTROFORCES_EIGEN_TARGET "")
if(TARGET Eigen3::Eigen)
  set(ASTROFORCES_EIGEN_TARGET Eigen3::Eigen)
elseif(TARGET eigen)
  set(ASTROFORCES_EIGEN_TARGET eigen)
endif()

add_library(astroforces_core INTERFACE)
target_include_directories(astroforces_core INTERFACE ${CMAKE_SOURCE_DIR}/libs/core/include)
target_link_libraries(astroforces_core INTERFACE fmt::fmt spdlog::spdlog)

add_library(astroforces_space_weather STATIC
  libs/space-weather/src/static_provider.cpp
  libs/space-weather/src/celestrak_csv_provider.cpp)
target_include_directories(astroforces_space_weather PUBLIC ${CMAKE_SOURCE_DIR}/libs/space-weather/include)
target_link_libraries(astroforces_space_weather PUBLIC astroforces_core)

add_library(astroforces_sc_props STATIC
  libs/sc-props/src/spacecraft.cpp)
target_include_directories(astroforces_sc_props PUBLIC ${CMAKE_SOURCE_DIR}/libs/sc-props/include)
target_link_libraries(astroforces_sc_props PUBLIC astroforces_core)

add_library(astroforces_models_basic STATIC
  libs/models-basic/src/exponential_atmosphere.cpp
  libs/models-basic/src/zero_wind.cpp)
target_include_directories(astroforces_models_basic PUBLIC ${CMAKE_SOURCE_DIR}/libs/models-basic/include)
target_link_libraries(astroforces_models_basic PUBLIC astroforces_core)

add_library(astroforces_forces_core STATIC
  libs/forces-core/src/core/perturbation.cpp
  libs/forces-core/src/gravity/gravity_sph_model.cpp
  libs/forces-core/src/gravity/tides/solid_earth_tide.cpp
  libs/forces-core/src/gravity/tides/pole_tide.cpp
  libs/forces-core/src/gravity/tides/constituent_tide.cpp
  libs/forces-core/src/gravity/tides/aod1b_tide.cpp
  libs/forces-core/src/gravity/tides/solid_earth_tide_freqdep.cpp
  libs/forces-core/src/gravity/gravity_sph_perturbation.cpp
  libs/forces-core/src/gravity/relativity_model.cpp
  libs/forces-core/src/gravity/relativity_perturbation.cpp
  libs/forces-core/src/gravity/third_body.cpp
  libs/forces-core/src/surface/drag/drag_model.cpp
  libs/forces-core/src/surface/drag/drag_perturbation.cpp
  libs/forces-core/src/surface/antenna_thrust/antenna_thrust_model.cpp
  libs/forces-core/src/surface/antenna_thrust/antenna_thrust_perturbation.cpp
  libs/forces-core/src/surface/earth_radiation/earth_radiation_model.cpp
  libs/forces-core/src/surface/earth_radiation/earth_radiation_perturbation.cpp
  libs/forces-core/src/surface/srp/srp_model.cpp
  libs/forces-core/src/surface/srp/srp_perturbation.cpp)
target_include_directories(astroforces_forces_core PUBLIC ${CMAKE_SOURCE_DIR}/libs/forces-core/include)
target_link_libraries(astroforces_forces_core PUBLIC astroforces_core astroforces_sc_props)
if(ASTROFORCES_EIGEN_TARGET)
  target_link_libraries(astroforces_forces_core PRIVATE ${ASTROFORCES_EIGEN_TARGET})
endif()

add_library(astroforces_adapters STATIC
  libs/adapters/src/nrlmsis21_adapter.cpp
  libs/adapters/src/dtm2020_adapter.cpp
  libs/adapters/src/hwm14_adapter.cpp)
target_include_directories(astroforces_adapters PUBLIC ${CMAKE_SOURCE_DIR}/libs/adapters/include)
target_link_libraries(astroforces_adapters PUBLIC
  astroforces_core
  msis21::msis21
  dtm2020::dtm2020
  hwm14::hwm14)

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

target_link_libraries(astroforces_forces_core PRIVATE jpl::eph)

add_executable(drag_cli apps/forces-cli/drag_cli.cpp)
target_link_libraries(drag_cli PRIVATE
  astroforces_core
  astroforces_space_weather
  astroforces_adapters
  astroforces_models_basic
  astroforces_forces_core)

add_executable(drag_batch_cli apps/forces-cli/drag_batch_cli.cpp)
target_link_libraries(drag_batch_cli PRIVATE
  astroforces_core
  astroforces_space_weather
  astroforces_adapters
  astroforces_models_basic
  astroforces_forces_core)

add_executable(perturbation_profile_cli apps/forces-cli/perturbation_profile_cli.cpp)
target_link_libraries(perturbation_profile_cli PRIVATE
  astroforces_core
  astroforces_space_weather
  astroforces_adapters
  astroforces_models_basic
  astroforces_forces_core)

add_executable(third_body_cli apps/forces-cli/third_body_cli.cpp)
target_link_libraries(third_body_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(third_body_batch_cli apps/forces-cli/third_body_batch_cli.cpp)
target_link_libraries(third_body_batch_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(srp_cli apps/forces-cli/srp_cli.cpp)
target_link_libraries(srp_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(srp_batch_cli apps/forces-cli/srp_batch_cli.cpp)
target_link_libraries(srp_batch_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(earth_radiation_cli apps/forces-cli/earth_radiation_cli.cpp)
target_link_libraries(earth_radiation_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(earth_radiation_batch_cli apps/forces-cli/earth_radiation_batch_cli.cpp)
target_link_libraries(earth_radiation_batch_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(antenna_thrust_cli apps/forces-cli/antenna_thrust_cli.cpp)
target_link_libraries(antenna_thrust_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(antenna_thrust_batch_cli apps/forces-cli/antenna_thrust_batch_cli.cpp)
target_link_libraries(antenna_thrust_batch_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(relativity_cli apps/forces-cli/relativity_cli.cpp)
target_link_libraries(relativity_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(gravity_sph_cli apps/forces-cli/gravity_sph_cli.cpp)
target_link_libraries(gravity_sph_cli PRIVATE
  astroforces_core
  astroforces_forces_core)

add_executable(frame_transform_cli apps/forces-cli/frame_transform_cli.cpp)
target_link_libraries(frame_transform_cli PRIVATE
  astroforces_core)

add_executable(transform_profile_cli apps/forces-cli/transform_profile_cli.cpp)
target_link_libraries(transform_profile_cli PRIVATE
  astroforces_core)

foreach(tgt astroforces_space_weather astroforces_sc_props astroforces_models_basic astroforces_forces_core astroforces_adapters drag_cli drag_batch_cli perturbation_profile_cli third_body_cli third_body_batch_cli srp_cli srp_batch_cli earth_radiation_cli earth_radiation_batch_cli antenna_thrust_cli antenna_thrust_batch_cli relativity_cli gravity_sph_cli frame_transform_cli transform_profile_cli)
  if(MSVC)
    target_compile_options(${tgt} PRIVATE /W4)
    if(ASTROFORCES_WARNINGS_AS_ERRORS)
      target_compile_options(${tgt} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${tgt} PRIVATE -Wall -Wextra -Wpedantic)
    if(ASTROFORCES_WARNINGS_AS_ERRORS)
      target_compile_options(${tgt} PRIVATE -Werror)
    endif()
  endif()
endforeach()

if(ASTROFORCES_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_executable(astroforces_test_drag_core tests/test_drag_core.cpp)
  target_link_libraries(astroforces_test_drag_core PRIVATE
    astroforces_core
    astroforces_space_weather
    astroforces_models_basic
    astroforces_forces_core)
  add_test(NAME DragCore COMMAND astroforces_test_drag_core)

  add_executable(astroforces_test_end_to_end tests/test_end_to_end.cpp)
  target_link_libraries(astroforces_test_end_to_end PRIVATE
    astroforces_core
    astroforces_space_weather
    astroforces_models_basic
    astroforces_forces_core)
  add_test(NAME DragEndToEnd COMMAND astroforces_test_end_to_end)

  add_executable(astroforces_test_adapters tests/test_adapters_integration.cpp)
  target_link_libraries(astroforces_test_adapters PRIVATE
    astroforces_core
    astroforces_space_weather
    astroforces_forces_core
    astroforces_adapters)
  target_compile_definitions(astroforces_test_adapters PRIVATE
    ASTROFORCES_NRLMSIS21_SOURCE_DIR="${nrlmsis21_ext_SOURCE_DIR}"
    ASTROFORCES_DTM2020_SOURCE_DIR="${dtm2020_ext_SOURCE_DIR}"
    ASTROFORCES_HWM14_SOURCE_DIR="${hwm14_ext_SOURCE_DIR}")
  add_test(NAME DragAdaptersIntegration COMMAND astroforces_test_adapters)

  add_executable(astroforces_test_space_weather tests/test_space_weather_celestrak.cpp)
  target_link_libraries(astroforces_test_space_weather PRIVATE
    astroforces_core
    astroforces_space_weather)
  add_test(NAME SpaceWeatherCelesTrak COMMAND astroforces_test_space_weather)

  add_executable(astroforces_test_space_weather_real_snippet tests/test_space_weather_celestrak_real_snippet.cpp)
  target_link_libraries(astroforces_test_space_weather_real_snippet PRIVATE
    astroforces_core
    astroforces_space_weather)
  target_compile_definitions(astroforces_test_space_weather_real_snippet PRIVATE
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  add_test(NAME SpaceWeatherCelesTrakRealSnippet COMMAND astroforces_test_space_weather_real_snippet)

  add_executable(astroforces_test_perturbation_stack tests/test_perturbation_stack.cpp)
  target_link_libraries(astroforces_test_perturbation_stack PRIVATE
    astroforces_core
    astroforces_space_weather
    astroforces_models_basic
    astroforces_forces_core)
  add_test(NAME PerturbationStack COMMAND astroforces_test_perturbation_stack)

  add_executable(astroforces_perf_benchmark tests/drag_perf_benchmark.cpp)
  target_link_libraries(astroforces_perf_benchmark PRIVATE
    astroforces_core
    astroforces_space_weather
    astroforces_models_basic
    astroforces_forces_core)

  add_executable(astroforces_test_third_body tests/test_third_body.cpp)
  target_link_libraries(astroforces_test_third_body PRIVATE
    astroforces_core
    astroforces_forces_core)
  target_compile_definitions(astroforces_test_third_body PRIVATE
    ASTROFORCES_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}"
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    $<$<BOOL:${ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA}>:ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA=1>)
  add_test(NAME ThirdBody COMMAND astroforces_test_third_body)

  add_executable(astroforces_test_srp tests/test_srp.cpp)
  target_link_libraries(astroforces_test_srp PRIVATE
    astroforces_core
    astroforces_forces_core)
  target_compile_definitions(astroforces_test_srp PRIVATE
    ASTROFORCES_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}"
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    $<$<BOOL:${ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA}>:ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA=1>)
  add_test(NAME SRP COMMAND astroforces_test_srp)

  add_executable(astroforces_test_earth_radiation tests/test_earth_radiation.cpp)
  target_link_libraries(astroforces_test_earth_radiation PRIVATE
    astroforces_core
    astroforces_forces_core)
  add_test(NAME EarthRadiation COMMAND astroforces_test_earth_radiation)

  add_executable(astroforces_test_antenna_thrust tests/test_antenna_thrust.cpp)
  target_link_libraries(astroforces_test_antenna_thrust PRIVATE
    astroforces_core
    astroforces_forces_core)
  add_test(NAME AntennaThrust COMMAND astroforces_test_antenna_thrust)

  add_executable(astroforces_test_relativity tests/test_relativity.cpp)
  target_link_libraries(astroforces_test_relativity PRIVATE
    astroforces_core
    astroforces_forces_core)
  target_compile_definitions(astroforces_test_relativity PRIVATE
    ASTROFORCES_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}"
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    $<$<BOOL:${ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA}>:ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA=1>)
  add_test(NAME Relativity COMMAND astroforces_test_relativity)

  add_executable(astroforces_test_gravity_sph tests/test_gravity_sph.cpp)
  target_link_libraries(astroforces_test_gravity_sph PRIVATE
    astroforces_core
    astroforces_forces_core)
  target_compile_definitions(astroforces_test_gravity_sph PRIVATE
    ASTROFORCES_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}"
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}"
    $<$<BOOL:${ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA}>:ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA=1>)
  add_test(NAME GravitySph COMMAND astroforces_test_gravity_sph)

  add_executable(astroforces_test_core_transforms tests/test_core_transforms.cpp)
  target_link_libraries(astroforces_test_core_transforms PRIVATE
    astroforces_core)
  add_test(NAME CoreTransforms COMMAND astroforces_test_core_transforms)

  add_executable(astroforces_test_eop_reader tests/test_eop_reader.cpp)
  target_link_libraries(astroforces_test_eop_reader PRIVATE
    astroforces_core)
  target_compile_definitions(astroforces_test_eop_reader PRIVATE
    ASTROFORCES_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  add_test(NAME EopReader COMMAND astroforces_test_eop_reader)

  add_executable(astroforces_test_cip_reader tests/test_cip_reader.cpp)
  target_link_libraries(astroforces_test_cip_reader PRIVATE
    astroforces_core)
  add_test(NAME CipReader COMMAND astroforces_test_cip_reader)

  add_executable(astroforces_test_leap_seconds tests/test_leap_seconds.cpp)
  target_link_libraries(astroforces_test_leap_seconds PRIVATE
    astroforces_core)
  add_test(NAME LeapSeconds COMMAND astroforces_test_leap_seconds)
endif()
