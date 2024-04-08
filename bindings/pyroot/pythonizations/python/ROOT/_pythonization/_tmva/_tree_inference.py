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


def SaveXGBoost(xgb_model, key_name, output_path, num_inputs=None, tmp_path="/tmp"):
    # Extract objective
    objective_map = {
        "multi:softprob": "softmax",  # Naming the objective softmax is more common today
        "binary:logistic": "logistic",
        "reg:linear": "identity",
        "reg:squarederror": "identity",
    }
    model_objective = xgb_model.objective
    if not model_objective in objective_map:
        raise Exception(
            'XGBoost model has unsupported objective "{}". Supported objectives are {}.'.format(
                model_objective, objective_map.keys()
            )
        )
    objective = cppyy.gbl.std.string(objective_map[model_objective])

    # Extract max depth of the trees
    max_depth = xgb_model.max_depth

    # Determine number of outputs
    if "reg:" in model_objective:
        num_outputs = 1
    elif "binary:" in model_objective:
        num_outputs = 1
    else:
        num_outputs = xgb_model.n_classes_

    # Dump XGB model to the tmp folder as json file
    import os
    import uuid

    tmp_path = os.path.join(tmp_path, str(uuid.uuid4()) + ".json")
    xgb_model.get_booster().dump_model(tmp_path, dump_format="json")

    import json

    with open(tmp_path, "r") as json_file:
        forest = json.load(json_file)

    bias = 0.0

    # Determine whether the model has a bias paramter and write bias trees
    if hasattr(xgb_model, "base_score") and "reg:" in model_objective:
        bias = xgb_model.base_score
        if not bias == 0.0:
            forest += [{"leaf": bias}] * num_outputs
    # print(str(forest).replace("u'", "'").replace("'", '"'))

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
    bdt = cppyy.gbl.TMVA.Experimental.RBDT.load_txt(output_path, features, num_outputs)

    bdt.baseResponse_ = bias

    with cppyy.gbl.TFile.Open(output_path, "RECREATE") as tFile:
        tFile.WriteObject(bdt, key_name)
