#!/usr/bin/python

import pandas as pd
import re


def formula2dict(formula):
    formula_dict = {}

    Cp = re.compile("C([0-9]{0,3})").findall(formula)
    if len(Cp) != 1:
        print(formula)
        raise ValueError("No Carbon in this formula %s" % (formula))
    formula_dict["C"] = Cp[0]

    Hp = re.compile("H([0-9]{0,3})").findall(formula)
    if len(Hp) != 1:
        print(formula)
        raise ValueError("No Hydrogen in this formula %s" % (formula))
    formula_dict["H"] = Hp[0]

    Np = re.compile("N([0-9]{0,3})").findall(formula)
    if len(Np) == 1:
        formula_dict["N"] = Np[0]

    Pp = re.compile("P([0-9]{0,3})").findall(formula)
    if len(Pp) == 1:
        formula_dict["P"] = Pp[0]

    Op = re.compile("O([0-9]{0,3})").findall(formula)
    if len(Op) == 1:
        formula_dict["O"] = Op[0]

#    Hap = re.compile("[Cl|Br|I|F]([0-9]{0,3})").findall(formula)
#    if len(Op) == 1:
#        formula_dict["O"] = Op[0]

    return(formula_dict)

def convert_dicts_to_integers(formula_dict_list):
    for f in formula_dict_list:
        for k,v in f.items():
            if v == "":
                v = f[k] = "1"
            f[k] = int(v)
    return(formula_dict_list)

def calculate_cOX(formula_dict):
    points = {"H": 1, "N": -3, "O": -2, "P": 5, "S": -2}
    sum = 0
    for k, v in formula_dict.items():
        if k == "C":
            continue
        sum += points[k] * v
    return(-sum/formula_dict["C"])

data = pd.read_csv("SIRIUS_CSI_FingerID_structures_annotations.csv", sep = "\t", nrows = 3)

formula_list = list(data["molecularFormula"].values)

dicts_from_formulas = [formula2dict(formula) for formula in formula_list]
dicts_from_formulas = convert_dicts_to_integers(dicts_from_formulas)

for idx, formula in zip(data["id"].values, dicts_from_formulas):
    print(idx, calculate_cOX(formula))
