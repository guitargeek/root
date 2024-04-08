# Author: Stefan Wunsch CERN  09/2019

################################################################################
# Copyright (C) 1995-2019, Rene Brun and Fons Rademakers.                      #
# All rights reserved.                                                         #
#                                                                              #
# For the licensing terms see $ROOTSYS/LICENSE.                                #
# For the list of contributors see $ROOTSYS/README/CREDITS.                    #
################################################################################

from .. import pythonization
import cppyy


def SaveXGBoost(xgb_model, key_name, output_path, num_inputs=None):

    # Determine number of input variables
    if not num_inputs is None:
        pass
    elif hasattr(xgb_model, "_features_count"):
        num_inputs = xgb_model._features_count
    else:
        raise Exception(
            "Failed to get number of input variables from XGBoost model. Please provide the additional keyword argument 'num_inputs' to this function."
        )

    xgb_model._Booster.dump_model(output_path)

    features = cppyy.gbl.std.vector["std::string"]([f"f{i}" for i in range(num_inputs)])
    bdt = cppyy.gbl.fastforest.load_txt(output_path, features)

    with cppyy.gbl.TFile.Open(output_path, "RECREATE") as tFile:
        tFile.WriteObject(bdt, key_name)
