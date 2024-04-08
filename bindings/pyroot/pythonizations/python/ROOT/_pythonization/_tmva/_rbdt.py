# Author: Stefan Wunsch CERN  09/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from .. import pythonization
from cppyy import gbl as gbl_namespace

def Compute(self, x):
    import numpy as np

    out = np.zeros(len(x))
    for i in range(len(x)):
        v = gbl_namespace.std.vector["float"](x[i])
        out[i] = self(v.data(), 0.0)
    return 1.0 / (1.0 + np.exp(-out))


@pythonization("FastForest", ns="fastforest", is_prefix=True)
def pythonize_rbdt(klass):
    # Parameters:
    # klass: class to be pythonized

    klass.Compute = Compute
