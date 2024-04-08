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


def get_basescore(model):
    import json

    """Get base score from an XGBoost sklearn estimator.

    Copy-pasted from XGBoost unit test code.
    """
    base_score = float(json.loads(model.get_booster().save_config())["learner"]["learner_model_param"]["base_score"])
    return base_score


def SaveXGBoost(xgb_model, key_name, output_path, num_inputs=None):

    import json

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
    num_outputs = 1
    if "multi:" in model_objective:
        num_outputs = xgb_model.n_classes_

    # Dump XGB model as json file
    xgb_model.get_booster().dump_model(output_path, dump_format="json")

    with open(output_path, "r") as json_file:
        forest = json.load(json_file)

    # Determine number of input variables
    if num_inputs is None:
        raise Exception(
            "Failed to get number of input variables from XGBoost model. Please provide the additional keyword argument 'num_inputs' to this function."
        )

    xgb_model._Booster.dump_model(output_path)

    features = cppyy.gbl.std.vector["std::string"]([f"f{i}" for i in range(num_inputs)])
    bdt = cppyy.gbl.TMVA.Experimental.RBDT.load_txt(output_path, features, num_outputs)

    bdt.logistic_ = objective == "logistic"
    if not bdt.logistic_:
        bdt.baseScore_ = get_basescore(xgb_model)

    with cppyy.gbl.TFile.Open(output_path, "RECREATE") as tFile:
        tFile.WriteObject(bdt, key_name)
