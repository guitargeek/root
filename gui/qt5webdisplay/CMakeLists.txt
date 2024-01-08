# Copyright (C) 1995-2023, Rene Brun and Fons Rademakers.
# All rights reserved.
#
# For the licensing terms see $ROOTSYS/LICENSE.
# For the list of contributors see $ROOTSYS/README/CREDITS.

############################################################################
# CMakeLists.txt file for building ROOT gui/qt5webdisplay package
############################################################################

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

find_package(Qt5 COMPONENTS Core WebEngine WebEngineWidgets CONFIG)

if(NOT Qt5_FOUND)
  message(FATAL_ERROR "Could NOT find Qt5 (WebEngine, WebEngineWidgets)")
endif()

set(CMAKE_AUTOMOC ON)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

ROOT_LINKER_LIBRARY(ROOTQt5WebDisplay
  rootwebview.cpp
  rootwebpage.cpp
  rooturlschemehandler.cpp
  rootqt5.cpp
  LIBRARIES
    Qt5::WebEngine
    Qt5::WebEngineWidgets
  DEPENDENCIES
    Core
    RHTTP
    ROOTWebDisplay
)
