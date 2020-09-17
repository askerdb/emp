#!/usr/bin/python

import pandas as pd
import re


def formula2dict(formula):
    formula_dict = {}

    c = re.compile("(C(?!l)|Cl|H|N|S|P|O|F|Br|I)([0-9]{0,3})").findall(formula)
    for i in c:
        formula_dict[i[0]] = i[1]

    #Check if there are formulas with no carbon and find out how to handle that case
    if "C" not in formula_dict:
        raise ValueError("No Carbon found in formula")
    
    for k,v in formula_dict.items():
        if v == "":
            v = formula_dict[k] = "1"

        formula_dict[k] = int(v)
    
    return(formula_dict)


def calculate_cOX(formula_dict):
    points = {"H": 1, "N": -3, "O": -2, "P": 5, "S": -2, "Cl": -1, "I": -1, "Br": -1, "F": -1}
    sum = 0
    for k, v in formula_dict.items():
        if k == "C":
            continue
        sum += points[k] * v
    return(-sum/formula_dict["C"])

data = pd.read_csv("SIRIUS_CSI_FingerID_structures_annotations.csv", sep = "\t")

formula_list = list(data["molecularFormula"].values)

dicts_from_formulas = [formula2dict(formula) for formula in formula_list]


for idx, formula in zip(data["id"].values, dicts_from_formulas):
    print(idx, calculate_cOX(formula))
