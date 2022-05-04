#!/usr/bin/env python3
# Author : Joana L. Rocha
import numpy as np
import pandas as pd
import sys
import os
from scipy.stats import chi2_contingency

infile = sys.argv[1]
outfile=sys.argv[2]
window_size = int(sys.argv[3])
step_size = int(sys.argv[4])
#window_type=sys.argv[5]

#if window_type not in ['snp', 'bp']:
#   print("please specify if you want to do windows on snps or basepairs with option snp or bp")

#tmp_outfile = outfile + '.unsorted'

offsets = np.arange(0, window_size, step_size)

np.set_printoptions(precision=4)

data = pd.read_csv(infile, '\t', names=["chromo", "position", "step", "LLratio", "global", "local", "fpop0", "fpop1", "fpop2", "fpop3", "fpop4"])

res = []

for offset in offsets:
    print(f"Running offset {offset}...")
    data["winstart"] = offset + window_size * ((data.position - offset) // window_size)
    data["winend"] = data.winstart + window_size
    #res_offset = data[data.winstart >= 0].groupby(["chromo", "winstart", "winend"]).LLratio.apply(lambda grp: grp.nlargest(500).sum()).to_frame()
    res_offset = data[data.winstart >= 0].groupby(["chromo", "winstart", "winend"]).LLratio.apply(lambda grp: grp.nlargest(500).mean()).to_frame()
    res.append(res_offset)
res = pd.concat(res)
res.to_csv(outfile, "\t", header=False)
