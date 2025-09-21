# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

# This package is only needed for graf3d implementation details! (TGLIncludes.h)

add_library(RGlew STATIC src/gl.c src/glx.c)
target_include_directories(RGlew PUBLIC inc/)
