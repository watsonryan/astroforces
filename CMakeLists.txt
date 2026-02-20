# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(astrodynamics_forces_cpp VERSION 0.1.0 LANGUAGES CXX)

option(DRAGCPP_BUILD_TESTS "Build tests" ON)
option(DRAGCPP_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)

set(DRAGCPP_NRLMSIS21_REPO "https://github.com/watsonryan/nrlmsis-2_1.git" CACHE STRING "NRLMSIS 2.1 repo")
set(DRAGCPP_DTM2020_REPO "https://github.com/watsonryan/dtm2020.git" CACHE STRING "DTM2020 repo")
set(DRAGCPP_HWM14_REPO "https://github.com/watsonryan/hwm14.git" CACHE STRING "HWM14 repo")
set(DRAGCPP_JPL_EPH_REPO "https://github.com/watsonryan/jplEphem.git" CACHE STRING "jplEphem repo")
set(DRAGCPP_EXTERNAL_MODELS_TAG "main" CACHE STRING "Branch/tag for external model repos")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(${CMAKE_SOURCE_DIR}/cmake/CPM.cmake)
set(DRAGCPP_EIGEN_REPO "https://gitlab.com/libeigen/eigen.git" CACHE STRING "Eigen repository URL")
set(DRAGCPP_EIGEN_TAG "5.0.0" CACHE STRING "Eigen git tag")

CPMAddPackage(
  NAME eigen
  GIT_REPOSITORY ${DRAGCPP_EIGEN_REPO}
  GIT_TAG ${DRAGCPP_EIGEN_TAG}
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

set(DRAGCPP_EIGEN_TARGET "")
if(TARGET Eigen3::Eigen)
  set(DRAGCPP_EIGEN_TARGET Eigen3::Eigen)
elseif(TARGET eigen)
  set(DRAGCPP_EIGEN_TARGET eigen)
endif()

add_library(dragcpp_atmo_core INTERFACE)
target_include_directories(dragcpp_atmo_core INTERFACE ${CMAKE_SOURCE_DIR}/libs/atmo-core/include)
target_link_libraries(dragcpp_atmo_core INTERFACE fmt::fmt spdlog::spdlog)

add_library(dragcpp_space_weather STATIC
  libs/space-weather/src/static_provider.cpp
  libs/space-weather/src/celestrak_csv_provider.cpp)
target_include_directories(dragcpp_space_weather PUBLIC ${CMAKE_SOURCE_DIR}/libs/space-weather/include)
target_link_libraries(dragcpp_space_weather PUBLIC dragcpp_atmo_core)

add_library(dragcpp_sc_props STATIC
  libs/sc-props/src/spacecraft.cpp)
target_include_directories(dragcpp_sc_props PUBLIC ${CMAKE_SOURCE_DIR}/libs/sc-props/include)
target_link_libraries(dragcpp_sc_props PUBLIC dragcpp_atmo_core)

add_library(dragcpp_models_basic STATIC
  libs/models-basic/src/exponential_atmosphere.cpp
  libs/models-basic/src/zero_wind.cpp)
target_include_directories(dragcpp_models_basic PUBLIC ${CMAKE_SOURCE_DIR}/libs/models-basic/include)
target_link_libraries(dragcpp_models_basic PUBLIC dragcpp_atmo_core)

add_library(dragcpp_drag_core STATIC
  libs/forces-core/src/drag_model.cpp
  libs/forces-core/src/drag_perturbation.cpp
  libs/forces-core/src/perturbation.cpp
  libs/forces-core/src/third_body.cpp
  libs/forces-core/src/srp_model.cpp
  libs/forces-core/src/srp_perturbation.cpp)
target_include_directories(dragcpp_drag_core PUBLIC ${CMAKE_SOURCE_DIR}/libs/forces-core/include)
target_link_libraries(dragcpp_drag_core PUBLIC dragcpp_atmo_core dragcpp_sc_props)
if(DRAGCPP_EIGEN_TARGET)
  target_link_libraries(dragcpp_drag_core PRIVATE ${DRAGCPP_EIGEN_TARGET})
endif()

add_library(dragcpp_adapters STATIC
  libs/adapters/src/placeholders.cpp
  libs/adapters/src/nrlmsis21_adapter.cpp
  libs/adapters/src/dtm2020_adapter.cpp
  libs/adapters/src/hwm14_adapter.cpp)
target_include_directories(dragcpp_adapters PUBLIC ${CMAKE_SOURCE_DIR}/libs/adapters/include)
target_link_libraries(dragcpp_adapters PUBLIC
  dragcpp_atmo_core
  msis21::msis21
  dtm2020::dtm2020
  hwm14::hwm14)

CPMAddPackage(
  NAME nrlmsis21_ext
  GIT_REPOSITORY ${DRAGCPP_NRLMSIS21_REPO}
  GIT_TAG ${DRAGCPP_EXTERNAL_MODELS_TAG}
  OPTIONS
    "MSIS21_BUILD_TESTS OFF"
    "MSIS21_BUILD_CLI OFF"
)

CPMAddPackage(
  NAME dtm2020_ext
  GIT_REPOSITORY ${DRAGCPP_DTM2020_REPO}
  GIT_TAG ${DRAGCPP_EXTERNAL_MODELS_TAG}
  OPTIONS
    "DTM2020_BUILD_TESTS OFF"
    "DTM2020_ENABLE_RESEARCH OFF"
)

CPMAddPackage(
  NAME hwm14_ext
  GIT_REPOSITORY ${DRAGCPP_HWM14_REPO}
  GIT_TAG ${DRAGCPP_EXTERNAL_MODELS_TAG}
  OPTIONS
    "HWM14_BUILD_TESTS OFF"
    "HWM14_BUILD_EXAMPLES OFF"
)

CPMAddPackage(
  NAME jpl_eph_ext
  GIT_REPOSITORY ${DRAGCPP_JPL_EPH_REPO}
  GIT_TAG ${DRAGCPP_EXTERNAL_MODELS_TAG}
  OPTIONS
    "JPL_EPH_BUILD_TESTS OFF"
    "JPL_EPH_BUILD_TOOLS OFF"
    "JPL_EPH_FETCH_DEPS OFF"
)

target_link_libraries(dragcpp_drag_core PRIVATE jpl::eph)

add_executable(drag_cli apps/drag-cli/main.cpp)
target_link_libraries(drag_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_space_weather
  dragcpp_sc_props
  dragcpp_adapters
  dragcpp_models_basic
  dragcpp_drag_core)

add_executable(drag_batch_cli apps/drag-cli/drag_batch_cli.cpp)
target_link_libraries(drag_batch_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_space_weather
  dragcpp_sc_props
  dragcpp_adapters
  dragcpp_models_basic
  dragcpp_drag_core)

add_executable(perturbation_profile_cli apps/drag-cli/perturbation_profile_cli.cpp)
target_link_libraries(perturbation_profile_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_space_weather
  dragcpp_sc_props
  dragcpp_adapters
  dragcpp_models_basic
  dragcpp_drag_core)

add_executable(third_body_cli apps/forces-cli/third_body_cli.cpp)
target_link_libraries(third_body_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_drag_core)

add_executable(third_body_batch_cli apps/forces-cli/third_body_batch_cli.cpp)
target_link_libraries(third_body_batch_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_drag_core)

add_executable(srp_cli apps/forces-cli/srp_cli.cpp)
target_link_libraries(srp_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_sc_props
  dragcpp_drag_core)

add_executable(srp_batch_cli apps/forces-cli/srp_batch_cli.cpp)
target_link_libraries(srp_batch_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_sc_props
  dragcpp_drag_core)

foreach(tgt dragcpp_space_weather dragcpp_sc_props dragcpp_models_basic dragcpp_drag_core dragcpp_adapters drag_cli drag_batch_cli perturbation_profile_cli third_body_cli third_body_batch_cli srp_cli srp_batch_cli)
  if(MSVC)
    target_compile_options(${tgt} PRIVATE /W4)
    if(DRAGCPP_WARNINGS_AS_ERRORS)
      target_compile_options(${tgt} PRIVATE /WX)
    endif()
  else()
    target_compile_options(${tgt} PRIVATE -Wall -Wextra -Wpedantic)
    if(DRAGCPP_WARNINGS_AS_ERRORS)
      target_compile_options(${tgt} PRIVATE -Werror)
    endif()
  endif()
endforeach()

if(DRAGCPP_BUILD_TESTS)
  include(CTest)
  enable_testing()
  add_executable(dragcpp_test_drag_core tests/test_drag_core.cpp)
  target_link_libraries(dragcpp_test_drag_core PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather
    dragcpp_sc_props
    dragcpp_models_basic
    dragcpp_drag_core)
  add_test(NAME DragCore COMMAND dragcpp_test_drag_core)

  add_executable(dragcpp_test_end_to_end tests/test_end_to_end.cpp)
  target_link_libraries(dragcpp_test_end_to_end PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather
    dragcpp_sc_props
    dragcpp_models_basic
    dragcpp_drag_core)
  add_test(NAME DragEndToEnd COMMAND dragcpp_test_end_to_end)

  add_executable(dragcpp_test_adapters tests/test_adapters_integration.cpp)
  target_link_libraries(dragcpp_test_adapters PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather
    dragcpp_sc_props
    dragcpp_drag_core
    dragcpp_adapters)
  target_compile_definitions(dragcpp_test_adapters PRIVATE
    DRAGCPP_NRLMSIS21_SOURCE_DIR="${nrlmsis21_ext_SOURCE_DIR}"
    DRAGCPP_DTM2020_SOURCE_DIR="${dtm2020_ext_SOURCE_DIR}"
    DRAGCPP_HWM14_SOURCE_DIR="${hwm14_ext_SOURCE_DIR}")
  add_test(NAME DragAdaptersIntegration COMMAND dragcpp_test_adapters)

  add_executable(dragcpp_test_space_weather tests/test_space_weather_celestrak.cpp)
  target_link_libraries(dragcpp_test_space_weather PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather)
  add_test(NAME SpaceWeatherCelesTrak COMMAND dragcpp_test_space_weather)

  add_executable(dragcpp_test_space_weather_real_snippet tests/test_space_weather_celestrak_real_snippet.cpp)
  target_link_libraries(dragcpp_test_space_weather_real_snippet PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather)
  target_compile_definitions(dragcpp_test_space_weather_real_snippet PRIVATE DRAGCPP_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
  add_test(NAME SpaceWeatherCelesTrakRealSnippet COMMAND dragcpp_test_space_weather_real_snippet)

  add_executable(dragcpp_test_perturbation_stack tests/test_perturbation_stack.cpp)
  target_link_libraries(dragcpp_test_perturbation_stack PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather
    dragcpp_sc_props
    dragcpp_models_basic
    dragcpp_drag_core)
  add_test(NAME PerturbationStack COMMAND dragcpp_test_perturbation_stack)

  add_executable(dragcpp_perf_benchmark tests/drag_perf_benchmark.cpp)
  target_link_libraries(dragcpp_perf_benchmark PRIVATE
    dragcpp_atmo_core
    dragcpp_space_weather
    dragcpp_sc_props
    dragcpp_models_basic
    dragcpp_drag_core)

  add_executable(dragcpp_test_third_body tests/test_third_body.cpp)
  target_link_libraries(dragcpp_test_third_body PRIVATE
    dragcpp_atmo_core
    dragcpp_sc_props
    dragcpp_drag_core)
  target_compile_definitions(dragcpp_test_third_body PRIVATE DRAGCPP_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}")
  add_test(NAME ThirdBody COMMAND dragcpp_test_third_body)

  add_executable(dragcpp_test_srp tests/test_srp.cpp)
  target_link_libraries(dragcpp_test_srp PRIVATE
    dragcpp_atmo_core
    dragcpp_sc_props
    dragcpp_drag_core)
  target_compile_definitions(dragcpp_test_srp PRIVATE DRAGCPP_JPL_EPH_SOURCE_DIR="${jpl_eph_ext_SOURCE_DIR}")
  add_test(NAME SRP COMMAND dragcpp_test_srp)
endif()
