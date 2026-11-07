# Authors:
# * Jonas Rembser 05/2021
# * Harshal Shende 06/2021

################################################################################
# Copyright (C) 1995-2020, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################


class RooAbsCollection(object):
    __cpp_name__ = 'RooAbsCollection'

    def addClone(self, arg, silent=False):

        import ROOT

        clonedArg = self._addClone(arg, silent)
        # There are two overloads of RooAbsCollection::addClone():
        #
        #   - RooAbsArg *addClone(const RooAbsArg& var, bool silent=false);
        #   - void addClone(const RooAbsCollection& list, bool silent=false);
        #
        # In the case of the RooAbsArg overload, we need to tell Python that it
        # doesn't own the returned pointer. That's because the function name
        # contains "Clone", which makes cppyy guess that the returned pointer
        # points to a clone owned by the caller. In the case of the
        # RooAbsCollection input, the return value will be `None` and we don't
        # need to change any ownership flags (in fact, calling
        # SetOwnership(None, False) would cause a crash).
        if clonedArg is not None:
            ROOT.SetOwnership(clonedArg, False)

    def addOwned(self, arg, silent=False):

        import ROOT

        self._addOwned(arg, silent)
        ROOT.SetOwnership(arg, False)
