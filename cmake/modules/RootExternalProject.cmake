# Copyright (C) 1995-2026, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

#----------------------------------------------------------------------------
# Downloaded builtins: the two-mirror policy and ROOT_DOWNLOAD_DIR
#
# Every source archive that ROOT downloads for a builtin_* option should be
# served by two independent locations that provide byte-identical files:
#
#   1. the CERN mirror, ${lcgpackages}/<archive> (or root.cern), and
#   2. an external mirror, normally the upstream release tarball on the
#      project's own site or the GitHub release asset or tag archive.
#
# Builtins that satisfy this are declared with ROOT_EXTERNALPROJECT_ADD below,
# listing both URLs. CMake tries the mirrors in order, so the external one also
# serves as a fallback when the CERN mirror is unavailable.
#
# Builtins that are not downloaded as an archive (git clones) or for which no
# external mirror serves the identical archive keep calling ExternalProject_Add
# directly, with a comment explaining why, and announce themselves with
# ROOT_REGISTER_NETWORK_ONLY_BUILTIN. They cannot be side-loaded through
# ROOT_DOWNLOAD_DIR.
#
# The archives are downloaded into ROOT_DOWNLOAD_DIR (see RootBuildOptions.cmake),
# named as on the CERN mirror. Archives found there with the expected hash are
# used instead of being downloaded, which allows to build ROOT without network
# access. ROOT_CHECK_BUILTIN_DOWNLOADS runs at the end of the configuration: it
# verifies the archives in ROOT_DOWNLOAD_DIR, and, if some are missing and there
# is no internet connection, stops with the list of archives that the enabled
# builtins need and where to get them.
#----------------------------------------------------------------------------

include(ExternalProject)

#----------------------------------------------------------------------------
# function ROOT_EXTERNALPROJECT_ADD(<name>
#                                   DOWNLOAD_DIR ${ROOT_DOWNLOAD_DIR}
#                                   URL <cern-mirror-url> <external-mirror-url>...
#                                   URL_HASH <algo>=<hash>
#                                   [<other ExternalProject_Add arguments>...])
#
# Drop-in replacement for ExternalProject_Add for builtins that follow the
# two-mirror policy described above. In addition to ExternalProject_Add it
#   - requires DOWNLOAD_DIR ${ROOT_DOWNLOAD_DIR}, at least two URLs and a
#     URL_HASH, and fails otherwise;
#   - names the downloaded archive after the first (CERN mirror) URL;
#   - registers the target for ROOT_CHECK_BUILTIN_DOWNLOADS.
#----------------------------------------------------------------------------
function(ROOT_EXTERNALPROJECT_ADD name)
  # Quoted, to keep empty arguments such as CONFIGURE_COMMAND "" (ExternalProject_Add
  # iterates over its arguments as a list, so it still sees them separately).
  set(_args "${ARGN}")

  set(_download_dir "")
  list(FIND _args DOWNLOAD_DIR _dir_index)
  if(_dir_index GREATER_EQUAL 0)
    math(EXPR _dir_index "${_dir_index} + 1")
    list(GET _args ${_dir_index} _download_dir)
  endif()
  if(NOT _download_dir STREQUAL ROOT_DOWNLOAD_DIR)
    message(FATAL_ERROR "ROOT_EXTERNALPROJECT_ADD(${name}): needs DOWNLOAD_DIR \${ROOT_DOWNLOAD_DIR}, "
      "so that the archive can be pre-downloaded there.")
  endif()

  list(FIND _args URL _url_index)
  if(_url_index LESS 0)
    message(FATAL_ERROR "ROOT_EXTERNALPROJECT_ADD(${name}): no URL given. "
      "Builtins that are not downloaded as an archive must use ExternalProject_Add directly.")
  endif()
  if(NOT "DOWNLOAD_NAME" IN_LIST _args)
    # Pin the local file name to the one of the CERN mirror, so that it does not
    # depend on the order of the mirrors or on ExternalProject's own rules.
    math(EXPR _url_index "${_url_index} + 1")
    list(GET _args ${_url_index} _first_url)
    get_filename_component(_download_name "${_first_url}" NAME)
    list(APPEND _args DOWNLOAD_NAME "${_download_name}")
  endif()

  ExternalProject_Add(${name} "${_args}")

  ExternalProject_Get_Property(${name} URL)
  list(LENGTH URL _nurls)
  if(_nurls LESS 2)
    message(FATAL_ERROR "ROOT_EXTERNALPROJECT_ADD(${name}): needs at least two URLs serving the identical "
      "archive, the CERN mirror first and then an external mirror, but got only: ${URL}\n"
      "Builtins without an external mirror must use ExternalProject_Add directly, "
      "see cmake/modules/RootExternalProject.cmake.")
  endif()
  # ExternalProject_Get_Property() aborts with a generic error when a property is
  # not set, so test for it first (the _EP_ prefix is what that function uses too).
  get_property(_has_hash TARGET ${name} PROPERTY _EP_URL_HASH SET)
  if(NOT _has_hash)
    message(FATAL_ERROR "ROOT_EXTERNALPROJECT_ADD(${name}): needs URL_HASH <algo>=<hash>, the two mirrors "
      "are only interchangeable because the hash pins the archive content.")
  endif()

  set_property(GLOBAL APPEND PROPERTY ROOT_EXTERNALPROJECT_TARGETS ${name})
endfunction()

#----------------------------------------------------------------------------
# function ROOT_REGISTER_NETWORK_ONLY_BUILTIN(<option>)
#
# To be called by builtins that cannot follow the two-mirror policy yet and
# therefore use ExternalProject_Add directly. Their archive can only be
# downloaded during the build, which ROOT_CHECK_BUILTIN_DOWNLOADS reports.
#----------------------------------------------------------------------------
function(ROOT_REGISTER_NETWORK_ONLY_BUILTIN option)
  set_property(GLOBAL APPEND PROPERTY ROOT_NETWORK_ONLY_BUILTINS ${option})
endfunction()

# Appends a description of the archive of an ExternalProject target declared
# with ROOT_EXTERNALPROJECT_ADD to the variable <outvar> in the caller's scope:
# file name, expected hash and all mirrors, optionally with a status after the
# file name.
function(_root_describe_archive target status outvar)
  ExternalProject_Get_Property(${target} URL URL_HASH DOWNLOAD_NAME)
  string(REGEX REPLACE "=.*$" "" _algo "${URL_HASH}")
  string(REGEX REPLACE "^[^=]*=" "" _expected "${URL_HASH}")
  string(TOLOWER "${_expected}" _expected)
  list(GET URL 0 _cern_mirror)
  list(SUBLIST URL 1 -1 _external_mirrors)
  set(_text "  ${DOWNLOAD_NAME} (${target})${status}\n")
  string(APPEND _text "    expected ${_algo}: ${_expected}\n")
  string(APPEND _text "    CERN mirror:     ${_cern_mirror}\n")
  foreach(_mirror IN LISTS _external_mirrors)
    string(APPEND _text "    external mirror: ${_mirror}\n")
  endforeach()
  set(${outvar} "${${outvar}}${_text}" PARENT_SCOPE)
endfunction()

#----------------------------------------------------------------------------
# function ROOT_CHECK_BUILTIN_DOWNLOADS()
#
# Called automatically at the end of the top-level configuration.
#
# Checks which archives of the targets declared with ROOT_EXTERNALPROJECT_ADD
# are already in ROOT_DOWNLOAD_DIR with the expected hash. The others, as well
# as the builtins registered with ROOT_REGISTER_NETWORK_ONLY_BUILTIN, have to be
# downloaded during the build. If the internet connectivity check did not run
# yet, it runs now to find out whether that is possible.
#
# Without internet connection, fails with a report of each missing or
# mismatching archive and how to provide it through ROOT_DOWNLOAD_DIR, naming
# also the enabled builtins that cannot be pre-downloaded and the options that
# were disabled for lack of connection.
#
# This failure is always fatal; the fail-on-missing option is kept for
# backward compatibility only and new code does not consult it.
#----------------------------------------------------------------------------
function(ROOT_CHECK_BUILTIN_DOWNLOADS)
  get_property(_targets GLOBAL PROPERTY ROOT_EXTERNALPROJECT_TARGETS)
  get_property(_network_only GLOBAL PROPERTY ROOT_NETWORK_ONLY_BUILTINS)
  get_property(_disabled GLOBAL PROPERTY ROOT_OPTIONS_DISABLED_OFFLINE)

  set(_verified)
  set(_to_download)
  set(_problems "")
  foreach(_target IN LISTS _targets)
    ExternalProject_Get_Property(${_target} URL_HASH DOWNLOAD_DIR DOWNLOAD_NAME)
    string(REGEX REPLACE "=.*$" "" _algo "${URL_HASH}")
    string(REGEX REPLACE "^[^=]*=" "" _expected "${URL_HASH}")
    string(TOLOWER "${_expected}" _expected)
    set(_file "${DOWNLOAD_DIR}/${DOWNLOAD_NAME}")
    if(NOT EXISTS "${_file}")
      _root_describe_archive(${_target} ": missing" _problems)
      list(APPEND _to_download ${DOWNLOAD_NAME})
    else()
      file(${_algo} "${_file}" _actual)
      if(_actual STREQUAL _expected)
        list(APPEND _verified ${DOWNLOAD_NAME})
      else()
        # ExternalProject replaces it during the build if it can download
        _root_describe_archive(${_target} ": ${_algo} mismatch, found ${_actual}" _problems)
        list(APPEND _to_download ${DOWNLOAD_NAME})
      endif()
    endif()
  endforeach()

  if(_verified)
    list(JOIN _verified " " _verified)
    message(STATUS "Archives verified in ${ROOT_DOWNLOAD_DIR}: ${_verified}")
  endif()

  # The connectivity check did not run if no option asked for it, or if it was
  # skipped because ROOT_DOWNLOAD_DIR already contains archives.
  get_property(_checked GLOBAL PROPERTY ROOT_NO_CONNECTION SET)
  if((_to_download OR _network_only) AND NOT _checked AND check_connection)
    ROOT_PROBE_CONNECTION()
  endif()
  get_property(_offline GLOBAL PROPERTY ROOT_NO_CONNECTION)

  if(NOT _offline)
    if(_to_download)
      list(JOIN _to_download " " _to_download)
      message(STATUS "Archives downloaded into ${ROOT_DOWNLOAD_DIR} during the build: ${_to_download}")
    endif()
    if(_network_only)
      list(JOIN _network_only " " _network_only)
      message(STATUS "Builtins without an external mirror, downloaded during the build regardless of ROOT_DOWNLOAD_DIR: ${_network_only}")
    endif()
    return()
  endif()

  set(_message "")
  if(_problems)
    string(APPEND _message
      "No internet connection, but ROOT_DOWNLOAD_DIR=${ROOT_DOWNLOAD_DIR} does not provide all archives needed by "
      "the enabled builtins:\n"
      "${_problems}"
      "Download each archive from either mirror (they serve identical files) and store it in ROOT_DOWNLOAD_DIR "
      "under the file name shown above, which is the one of the CERN mirror. You can also point "
      "-DROOT_DOWNLOAD_DIR=<dir> to another directory with these archives.\n")
  endif()
  if(_network_only)
    list(JOIN _network_only " " _network_only)
    string(APPEND _message
      "No internet connection, but these enabled builtins must download their archive during the build, "
      "because no external mirror provides it and they therefore cannot be pre-downloaded into ROOT_DOWNLOAD_DIR yet: "
      "${_network_only}. Disable them or use the corresponding system packages.\n")
  endif()
  if(_message)
    if(_disabled)
      list(JOIN _disabled " " _disabled)
      string(APPEND _message
        "Furthermore, the archives of options that were disabled for lack of connection are not listed above: "
        "${_disabled}. To have them listed too, re-enable them, together with -Dfail-on-missing=ON or after "
        "placing the archives above in ROOT_DOWNLOAD_DIR (which skips the connectivity check).\n")
    endif()
    message(FATAL_ERROR
      "${_message}"
      "Alternatively, check your connection, disable the builtins, or bypass the connection check with -Dcheck_connection=OFF.")
  elseif(_disabled)
    list(JOIN _disabled " " _disabled)
    message(STATUS "Options disabled for lack of internet connection: ${_disabled}. To build them without "
      "internet connection, re-enable them with -Dfail-on-missing=ON to get the list of their archives.")
  elseif(NOT _targets)
    # The connectivity check only runs for options that need the network, so
    # something needs it but was not registered: do not let the build fail late.
    message(FATAL_ERROR "No internet connection, but an enabled option requires one. Check your connection, "
      "or bypass the connection check with -Dcheck_connection=OFF.")
  endif()
endfunction()

cmake_language(DEFER DIRECTORY ${CMAKE_SOURCE_DIR} CALL ROOT_CHECK_BUILTIN_DOWNLOADS)
