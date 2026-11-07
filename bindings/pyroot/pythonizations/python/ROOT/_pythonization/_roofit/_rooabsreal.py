# Authors:
# * Hinnerk C. Schmidt 02/2021
# * Jonas Rembser 03/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


class RooAbsReal(object):
    __cpp_name__ = "RooAbsReal"

    def setEvalErrorLoggingMode(m):

        import ROOT

        if isinstance(m, str):
            # Hardcode enum integer values here, because enum lookups cause
            # some memory fiasco at the end.
            # TODO: fix this in cppyy / PyROOT!
            lut = {"PrintErrors": 0, "CollectErrors": 1, "CountErrors": 2, "Ignore": 3}
            if m not in lut:
                raise ValueError(
                    "Unsupported value passed. The value has to be the name or enum value of an RooAbsReal ErrorLogging mode: {}".format(
                        lut
                    )
                )

            m = lut[m]

        ROOT.RooAbsReal._setEvalErrorLoggingMode(m)
