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


r"""
/**
\class RooDataSet
\brief \parblock \endparblock
\htmlonly
<div class="pyrootbox">
\endhtmlonly

## PyROOT

Some member functions of RooDataSet that take a RooCmdArg as argument also support keyword arguments.
So far, this applies to RooDataSet() constructor and RooDataSet::plotOnXY.
For example, the following code is equivalent in PyROOT:
\code{.py}
# Directly passing a RooCmdArg:
dxy = ROOT.RooDataSet("dxy", "dxy", ROOT.RooArgSet(x, y), ROOT.RooFit.StoreError(ROOT.RooArgSet(x, y)))

# With keyword arguments:
dxy = ROOT.RooDataSet("dxy", "dxy", ROOT.RooArgSet(x, y), StoreError=(ROOT.RooArgSet(x, y)))
\endcode

\htmlonly
</div>
\endhtmlonly
*/
"""

from ._utils import _kwargs_to_roocmdargs


class RooDataSet(object):
    def __init__(self, *args, **kwargs):
        # Redefinition of `RooDataSet` constructor for keyword arguments.
        # The keywords must correspond to the CmdArg of the constructor function.
        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        self._init(*args, **kwargs)

    def plotOnXY(self, *args, **kwargs):
        # Redefinition of `RooDataSet.plotOnXY` for keyword arguments.
        # The keywords must correspond to the CmdArg of the `plotOnXY` function.
        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        return self._plotOnXY(*args, **kwargs)

    @staticmethod
    def from_numpy(data, name=None, title=None, variables=None):

        import ROOT

        import ctypes
        import numpy as np

        name = "" if name is None else name
        title = "" if title is None else title
        n_entries = None

        data_real = ROOT.std.vector["double*"]()

        for arr_name, arr in data.items():
            if arr.dtype != np.float64:
                raise TypeError("Arrays passed to RooDataSet.from_numpy() need to be of type float64.")
            if n_entries is None:
                n_entries = len(arr)
            elif n_entries != len(arr):
                raise RuntimeError("Arrays passed to RooDataSet.from_numpy() must be all of same length.")
            if len(arr.shape) != 1:
                raise RuntimeError("Arrays passed to RooDataSet.from_numpy() must be one dimensional.")

            arr_ptr = arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            data_real.push_back(arr_ptr)

        if variables is None:
            variables = []
            for arr_name, arr in data.items():
                variables.append(ROOT.RooRealVar(arr_name, arr_name, arr[0]))

        data_set = ROOT.RooDataSet.fromArrays(name, title, variables, n_entries, data_real)

        return data_set

    @staticmethod
    def from_pandas(df, name=None, title=None, variables=None):

        import ROOT

        data = {}
        for column in df:
            data[column] = df[column].values
        return ROOT.RooDataSet.from_numpy(data, name=name, title=title, variables=variables)
