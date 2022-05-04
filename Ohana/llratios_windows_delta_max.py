#!/usr/bin/env python3
# Author : Joana L. Rocha
import pandas as pd
import numpy as np
import sys

branch1= sys.argv[1]
branch2= sys.argv[2]
output= sys.argv[3]

df1=pd.read_csv(branch1, '\t', header=None,  names=["chromo", "pos_start", "pos_end", "llratio1", "nsites1"])

df2=pd.read_csv(branch2, '\t', header=None,  names=["chromo",  "pos_start", "pos_end", "llratio2", "nsites2"])

print(df1.head())
print(df2.head())

merged_llratios = df1.merge(df2, how = 'inner', on=["chromo", "pos_start", "pos_end"])
merged_llratios['max'] = merged_llratios[["llratio1", "llratio2"]].max(axis=1)
merged_llratios['delta'] = merged_llratios['llratio1'] - merged_llratios['llratio2']

print(merged_llratios.head())

merged_llratios.to_csv(output, sep='\t', header=False, index=False)
