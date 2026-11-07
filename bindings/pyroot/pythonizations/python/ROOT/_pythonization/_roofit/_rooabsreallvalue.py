# Authors:
# * Jonas Rembser 06/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


class RooAbsRealLValue(object):

    __cpp_name__ = 'RooAbsRealLValue'

    def frame(self, *args, **kwargs):
        frame = self._frame(*args, **kwargs)

        # Add a Python reference to the plot variable to the RooPlot because
        # the plot variable needs to survive at least as long as the plot.
        frame._plotVar = self

        return frame
