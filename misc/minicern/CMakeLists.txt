# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT misc/minicern package
############################################################################

ROOT_LINKER_LIBRARY(minicern *.c *.f TYPE STATIC)
set_property(TARGET minicern PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(minicern ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})

# Disable warnings, since what has worked for 40 years...
# (see https://sft.its.cern.ch/jira/browse/ROOT-9179 for the warnings)
#
# All remaining warnings are legacy F77 argument rank mismatches, where a scalar
# (often a literal like 0) is passed to a routine that declares the dummy as an
# array but only ever touches element (1) -- harmless storage association.
# gfortran offers no granular way to silence just this class (-Wno-argument-mismatch
# has no effect), so the blanket -w stays. The related "actual argument contains too
# few elements" warnings *were* fixed by giving those dummy arguments assumed-size
# (*) declarations, which is behaviour-preserving for the element-wise access here.
target_compile_options(minicern PRIVATE "$<$<COMPILE_LANGUAGE:Fortran>:-w>")

# The offset computations between common blocks go wrong an arm64 with -O2.
# This is not a proper fix, but getting it right would require rewriting zebra's
# distance between common blocks logic ...
target_compile_options(minicern PUBLIC "$<$<COMPILE_LANGUAGE:Fortran>:-O0>")