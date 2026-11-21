# Authors:
# * Jonas Rembser 11/2026

################################################################################
# Copyright (C) 1995-2026, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from .. import pythonization
from ._utils import _kwargs_to_roocmdargs, cpp_signature


class ModelConfig(object):
    r"""The RooFit::ModelConfig::fitTo() and RooFit::ModelConfig::createNLL() functions support keyword arguments,
    like the corresponding functions of RooAbsPdf.
    For example, the following code is equivalent in PyROOT:
    \code{.py}
    # Directly passing a RooCmdArg:
    modelConfig.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1))

    # With keyword arguments:
    modelConfig.fitTo(data, Save=True, PrintLevel=-1)
    \endcode"""

    __cpp_name__ = "RooStats::ModelConfig"

    @cpp_signature("RooFit::ModelConfig::fitTo()")
    def _ModelConfig_fitTo(self, *args, **kwargs):
        r"""The RooFit::ModelConfig::fitTo() function is pythonized with the command argument pythonization.
        The keywords must correspond to the CmdArgs of the function.
        """
        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        return self._fitTo(*args, **kwargs)

    @cpp_signature("RooFit::ModelConfig::createNLL()")
    def _ModelConfig_createNLL(self, *args, **kwargs):
        r"""The RooFit::ModelConfig::createNLL() function is pythonized with the command argument pythonization.
        The keywords must correspond to the CmdArgs of the function.
        """
        args, kwargs = _kwargs_to_roocmdargs(*args, **kwargs)
        return self._createNLL(*args, **kwargs)

    @pythonization("ModelConfig", ns="RooStats")
    def pythonize_modelconfig(klass):
        klass._fitTo = klass.fitTo
        klass.fitTo = ModelConfig._ModelConfig_fitTo

        klass._createNLL = klass.createNLL
        klass.createNLL = ModelConfig._ModelConfig_createNLL
