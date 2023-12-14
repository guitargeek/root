# This header-only library contains common inlined math functions that can be
# used in RooBatchCompute, the code generation, and the traditional RooFit
# libraries.

set (EXTRA_DICT_OPTS)
if (runtime_cxxmodules AND WIN32)
  set (EXTRA_DICT_OPTS NO_CXXMODULE)
endif()

ROOT_STANDARD_LIBRARY_PACKAGE(RooFitMath
  HEADERS
    RooFitMath.h
  NO_SOURCES
  LINKDEF
    LinkDef.h
  DEPENDENCIES
    MathCore
  ${EXTRA_DICT_OPTS}
)
