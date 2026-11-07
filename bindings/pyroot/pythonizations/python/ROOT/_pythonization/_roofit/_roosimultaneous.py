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


from ._utils import _dict_to_flat_map, cpp_signature


class RooSimultaneous(object):

    __cpp_name__ = 'RooSimultaneous'

    @cpp_signature(
        "RooSimultaneous(const char *name, const char *title,"
        "                std::map<std::string,RooAbsPdf*> pdfMap, RooAbsCategoryLValue& inIndexCat) ;"
    )
    def __init__(self, *args):
        r"""The RooSimultaneous constructor that takes a map of category names
        to PDFs is accepting a Python dictionary in Python.
        """
        if len(args) >= 3 and isinstance(args[2], dict):
            args = list(args)
            args[2] = _dict_to_flat_map(args[2], {"std::string": "RooAbsPdf*"})
        self._init(*args)
