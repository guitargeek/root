# Regression test driver: rootcling must resolve the C++ standard library
# against the compiler that drives the build (handed to it via -compilerI),
# not against an unrelated compiler that merely happens to be first on $PATH.
#
# See https://github.com/root-project/root/issues (rootcling picking up the
# wrong, often too old, system compiler -- e.g. failing to find <concepts> --
# when the build compiler is not on $PATH, a common situation with relocatable
# /cvmfs toolchains).
#
# Strategy: shadow every common compiler name on $PATH with a fake compiler
# that advertises an *empty* C++ include directory (so the standard library
# headers cannot be found through it), while passing the real C++ standard
# library directories via -compilerI -- exactly as ROOT_GENERATE_DICTIONARY
# does with CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES. Dictionary generation must
# then succeed by using the -compilerI directories.
#
# Required -D arguments:
#   ROOTCLING     - path to the rootcling executable
#   HEADER        - header to generate a dictionary for
#   LINKDEF       - the LinkDef.h
#   CXX_INC_DIRS  - '|'-separated compiler implicit include directories
#   WORKDIR       - scratch directory (created/cleaned by this script)

if(NOT ROOTCLING OR NOT HEADER OR NOT LINKDEF OR NOT WORKDIR)
  message(FATAL_ERROR "compilerIPathFallback.cmake: missing required -D argument")
endif()

# Fresh scratch area.
file(REMOVE_RECURSE "${WORKDIR}")
file(MAKE_DIRECTORY "${WORKDIR}")

# A "c++"-named but empty include directory: it matches cling's `++` filter but
# contains none of the standard library headers.
set(badinc "${WORKDIR}/wrong-toolchain/include/c++/0")
file(MAKE_DIRECTORY "${badinc}")

# Fake compilers shadowing whichever name cling tries to invoke (CLING_CXX_RLTV
# is relative, so $PATH order decides which binary is used).
set(fakebin "${WORKDIR}/fakebin")
file(MAKE_DIRECTORY "${fakebin}")
foreach(name g++ c++ clang++ gcc cc clang)
  set(script "${fakebin}/${name}")
  # Emit just enough of `<compiler> -xc++ -E -v` for cling's header query.
  file(WRITE "${script}"
"#!/bin/sh
cat 1>&2 <<'EOF'
#include <...> search starts here:
 ${badinc}
End of search list.
EOF
exit 0
")
  execute_process(COMMAND chmod +x "${script}")
endforeach()

# Build the -compilerI arguments from the real C++ standard library directories,
# keeping only the C++ ones (mirroring what rootcling itself filters on).
string(REPLACE "|" ";" _inc_dirs "${CXX_INC_DIRS}")
set(compilerI_args "")
foreach(d IN LISTS _inc_dirs)
  if(d MATCHES "c\\+\\+")
    list(APPEND compilerI_args "-compilerI${d}")
  endif()
endforeach()
if(NOT compilerI_args)
  message(FATAL_ERROR "No C++ standard library directories among CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES; "
                      "cannot run this test on this configuration.")
endif()

# Put the fake compiler first on $PATH and make sure no user override masks the
# behaviour under test.
set(ENV{PATH} "${fakebin}:$ENV{PATH}")
unset(ENV{CLING_CPPSYSINCL})

set(dict "${WORKDIR}/G__compilerIProbe.cxx")

execute_process(
  COMMAND ${ROOTCLING} -v2 -f ${dict} ${compilerI_args} ${HEADER} ${LINKDEF}
  RESULT_VARIABLE _rc
  OUTPUT_VARIABLE _out
  ERROR_VARIABLE _err)

message(STATUS "rootcling stdout:\n${_out}")
message(STATUS "rootcling stderr:\n${_err}")

if(NOT _rc EQUAL 0)
  message(FATAL_ERROR
    "rootcling failed (exit ${_rc}): it resolved the C++ standard library against the bogus "
    "$PATH compiler instead of the -compilerI directories.")
endif()
if(NOT EXISTS "${dict}")
  message(FATAL_ERROR "rootcling reported success but did not produce ${dict}.")
endif()

message(STATUS "PASS: dictionary generated from the -compilerI directories despite a bogus $PATH compiler.")
file(REMOVE_RECURSE "${WORKDIR}")
