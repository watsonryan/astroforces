# Author: Watosn
cmake_minimum_required(VERSION 3.25)
project(drag_cpp VERSION 0.1.0 LANGUAGES CXX)

option(DRAGCPP_BUILD_TESTS "Build tests" ON)
option(DRAGCPP_WARNINGS_AS_ERRORS "Treat warnings as errors" ON)

set(DRAGCPP_NRLMSIS21_REPO "https://github.com/watsonryan/nrlmsis-2_1.git" CACHE STRING "NRLMSIS 2.1 repo")
set(DRAGCPP_DTM2020_REPO "https://github.com/watsonryan/dtm2020.git" CACHE STRING "DTM2020 repo")
set(DRAGCPP_HWM14_REPO "https://github.com/watsonryan/hwm14.git" CACHE STRING "HWM14 repo")
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

set(DRAGCPP_EIGEN_TARGET "")
if(TARGET Eigen3::Eigen)
  set(DRAGCPP_EIGEN_TARGET Eigen3::Eigen)
elseif(TARGET eigen)
  set(DRAGCPP_EIGEN_TARGET eigen)
endif()

add_library(dragcpp_atmo_core INTERFACE)
target_include_directories(dragcpp_atmo_core INTERFACE ${CMAKE_SOURCE_DIR}/libs/atmo-core/include)

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
  libs/drag-core/src/drag_model.cpp)
target_include_directories(dragcpp_drag_core PUBLIC ${CMAKE_SOURCE_DIR}/libs/drag-core/include)
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

add_executable(drag_cli apps/drag-cli/main.cpp)
target_link_libraries(drag_cli PRIVATE
  dragcpp_atmo_core
  dragcpp_space_weather
  dragcpp_sc_props
  dragcpp_adapters
  dragcpp_models_basic
  dragcpp_drag_core)

foreach(tgt dragcpp_space_weather dragcpp_sc_props dragcpp_models_basic dragcpp_drag_core dragcpp_adapters drag_cli)
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
endif()
