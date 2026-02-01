# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

ROOT_PYTHON_PACKAGE(cppyy
  SOURCES
  ROOT/_cppyy/_stdcpp_fix.py
  ROOT/_cppyy/__init__.py
  ROOT/_cppyy/_cpython_cppyy.py
  ROOT/_cppyy/_pypy_cppyy.py
  ROOT/_cppyy/_pythonization.py
  ROOT/_cppyy/_typemap.py
  ROOT/_cppyy/_version.py
  ROOT/_cppyy/interactive.py
  ROOT/_cppyy/ll.py
  ROOT/_cppyy/numba_ext.py
  ROOT/_cppyy/reflex.py
  ROOT/_cppyy/types.py
  cppyy/__init__.py
)
