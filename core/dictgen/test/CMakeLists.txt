# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# @author Danilo Piparo CERN

ROOT_ADD_GTEST(dictgen_base dictgen_base.cxx LIBRARIES Core)

# Regression test: rootcling must resolve the C++ standard library against the
# build compiler (passed via -compilerI), not against an unrelated compiler that
# happens to be first on $PATH. The driver fakes such a wrong $PATH compiler, so
# it relies on the POSIX shell and is skipped on Windows.
if(NOT MSVC)
  string(REPLACE ";" "|" _compilerI_inc_dirs "${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES}")
  ROOT_ADD_TEST(rootcling-compilerI-path-fallback
    COMMAND ${CMAKE_COMMAND}
            -DROOTCLING=$<TARGET_FILE:rootcling>
            -DHEADER=${CMAKE_CURRENT_SOURCE_DIR}/compilerICxxHeader.h
            -DLINKDEF=${CMAKE_CURRENT_SOURCE_DIR}/compilerICxxHeader_LinkDef.h
            -DCXX_INC_DIRS=${_compilerI_inc_dirs}
            -DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/rootcling-compilerI-path-fallback
            -P ${CMAKE_CURRENT_SOURCE_DIR}/compilerIPathFallback.cmake
    DEPENDS rootcling
    LABELS core dictgen)
endif()
