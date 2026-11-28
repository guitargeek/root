# windows fails with message
# PyCool_rflx.cxx(824,1): error C2026: string too big, trailing characters truncated

if(NOT MSVC OR win_broken_tests)
  ROOTTEST_ADD_TEST(PyCool
                    COMMAND ${ROOT_genreflex_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/PyCool.h
                            --interpreteronly -I ${CMAKE_CURRENT_SOURCE_DIR} -I ${CMAKE_CURRENT_SOURCE_DIR}/CoralBase
                            --select=${CMAKE_CURRENT_SOURCE_DIR}/PyCool_selection.xml
                            --output=${CMAKE_CURRENT_BINARY_DIR}/PyCool_rflx.cxx
                    FIXTURES_SETUP PyCool-reflex-dict)

  # Standalone project that builds a static library, only to ensure that the
  # generated dictionary can be compiled. It is built by a test via
  # "ctest --build-and-test", so it creates no target in this project and the
  # test can run in parallel with the others.
  set(_pycoollib_project ${CMAKE_CURRENT_BINARY_DIR}/PyCoolLib_project)
  file(CONFIGURE OUTPUT ${_pycoollib_project}/src/CMakeLists.txt @ONLY CONTENT [[
cmake_minimum_required(VERSION @CMAKE_MINIMUM_REQUIRED_VERSION@)
project(PyCoolLib LANGUAGES CXX)
set(CMAKE_CXX_STANDARD @CMAKE_CXX_STANDARD@)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(MSVC)
  # ROOT's headers need ROOT's compiler flags, e.g. -Zc:__cplusplus
  set(CMAKE_CXX_FLAGS "@CMAKE_CXX_FLAGS@")
endif()

add_library(PyCoolLib STATIC "@CMAKE_CURRENT_BINARY_DIR@/PyCool_rflx.cxx")
target_include_directories(PyCoolLib PRIVATE "@CMAKE_BINARY_DIR@/include" "@CMAKE_CURRENT_SOURCE_DIR@")
]])

  ROOTTEST_TARGETNAME_FROM_FILE(_pycoollib_test PyCoolLib)
  ROOT_ADD_TEST(${_pycoollib_test}-build
    BUILD PyCoolLib
    SOURCE_DIR ${_pycoollib_project}/src
    BINARY_DIR ${_pycoollib_project}/build
    FIXTURES_REQUIRED PyCool-reflex-dict)
endif()
