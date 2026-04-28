# Author: Watson
if(NOT CPM_VERSION)
  set(CPM_VERSION 0.40.2)
endif()

if(CPM_SOURCE_CACHE)
  set(CPM_DOWNLOAD_LOCATION "${CPM_SOURCE_CACHE}/cpm/CPM_${CPM_VERSION}.cmake")
elseif(DEFINED ENV{CPM_SOURCE_CACHE})
  set(CPM_DOWNLOAD_LOCATION "$ENV{CPM_SOURCE_CACHE}/cpm/CPM_${CPM_VERSION}.cmake")
else()
  set(CPM_DOWNLOAD_LOCATION "${CMAKE_BINARY_DIR}/cmake/CPM_${CPM_VERSION}.cmake")
endif()

if(NOT EXISTS ${CPM_DOWNLOAD_LOCATION})
  file(DOWNLOAD
    https://github.com/cpm-cmake/CPM.cmake/releases/download/v${CPM_VERSION}/CPM.cmake
    ${CPM_DOWNLOAD_LOCATION}
    STATUS cpm_download_status)
endif()

if(EXISTS ${CPM_DOWNLOAD_LOCATION})
  include(${CPM_DOWNLOAD_LOCATION} OPTIONAL RESULT_VARIABLE cpm_include_result)
endif()

if(NOT COMMAND CPMAddPackage)
  include(FetchContent)
  function(CPMAddPackage)
    set(options)
    set(oneValueArgs NAME GITHUB_REPOSITORY VERSION)
    set(multiValueArgs OPTIONS)
    cmake_parse_arguments(CPM "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
    if(NOT CPM_NAME)
      message(FATAL_ERROR "CPMAddPackage fallback requires NAME")
    endif()
    if(NOT CPM_GITHUB_REPOSITORY)
      message(FATAL_ERROR "CPMAddPackage fallback requires GITHUB_REPOSITORY")
    endif()
    if(NOT CPM_VERSION)
      message(FATAL_ERROR "CPMAddPackage fallback requires VERSION")
    endif()

    foreach(opt IN LISTS CPM_OPTIONS)
      if(opt MATCHES "^([^ ]+) (.+)$")
        set(${CMAKE_MATCH_1} ${CMAKE_MATCH_2} CACHE INTERNAL "")
      endif()
    endforeach()

    set(cpm_git_tag "v${CPM_VERSION}")
    if(CPM_GITHUB_REPOSITORY STREQUAL "fmtlib/fmt")
      set(cpm_git_tag "${CPM_VERSION}")
    endif()
    FetchContent_Declare(
      ${CPM_NAME}
      GIT_REPOSITORY https://github.com/${CPM_GITHUB_REPOSITORY}.git
      GIT_TAG ${cpm_git_tag}
    )
    FetchContent_MakeAvailable(${CPM_NAME})
  endfunction()
endif()
