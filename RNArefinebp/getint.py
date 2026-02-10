#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 22:11:43 2024

@author: tsukiyamashou
"""

import pandas as pd

#path="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3D/Examples/expnew/RNAref3D_bp.csv"
#outpath="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3D/Examples/expnew/RNAref3D_bp.txt"

path="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3D/Examples/expnew/RNAref3D_stack.csv"
outpath="/Users/sho/Documents/project/RNArefiner/RNArefiner/stacking_prediction/RNAref3D/Examples/expnew/RNAref3D_stack.txt"

df = pd.read_csv(path)

df.loc[:, "base_i"] = df.loc[:, "base_i"].astype(str)
df.loc[:, "base_j"] = df.loc[:, "base_j"].astype(str)
df.loc[:, "score"] = df.loc[:, "score"].astype(str)
interact = (df.loc[:, ["base_i", "base_j", "score"]]).values.tolist()

interact = [",".join(interact[i]) for i in range(len(interact))]

with open(outpath, "w") as f:
    f.write("\n".join(interact))
































