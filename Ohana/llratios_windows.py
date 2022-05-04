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
window_type=sys.argv[5]

if window_type not in ['snp', 'bp']:
    print("please specify if you want to do windows on snps or basepairs with option snp or bp")

tmp_outfile = outfile + '.unsorted'

offsets = np.arange(0, window_size, step_size)

np.set_printoptions(precision=4)

data = pd.read_csv(infile, '\t', names=["chromo", "position", "step", "LLratio", "global", "local", "fpop0", "fpop1", "fpop2", "fpop3", "fpop4"])


llratio = data.LLratio

n_snps = data.shape[0]

cum_llratio = np.cumsum(llratio)


def calculate_statistics(chr, pos_start, pos_end, i_start, i_end, f):
    if i_start > i_end: return

    nsites= i_end - i_start + 1

    llvalue = cum_llratio[i_end] - (cum_llratio[i_start-1] if i_start>0 else 0)

    f.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(
        chr, pos_start, pos_end,
        llvalue,
        nsites
    ))


with open(tmp_outfile,'w') as f:
    for offset in offsets:
        if window_type == 'bp':
            i = 0
            pos_start = offset

            while True:
                chr = data.chromo[i]
                while True:
                    if i == n_snps:
                        i_start = 'end'
                        break
                    elif data.chromo[i] != chr:
                        chr = data.chromo[i]
                        pos_start = offset
                    elif data.position[i] >= pos_start:
                        i_start = i
                        break
                    else:
                        i += 1

                if i_start == 'end': break

                i = i_start
                chr = data.chromo[i_start]
                pos_end = pos_start+window_size

                while True:
                    if i == n_snps:
                        i_end = 'end'
                        status = False
                        break
                    if data.chromo[i] != chr:
                        i_end = i
                        status = False
                        break
                    elif data.position[i] >= pos_end:
                        i_end = i
                        status = True
                        break
                    else:
                        i += 1

                if i_end == 'end': break
                if status:
                    calculate_statistics(data.chromo[i_start], pos_start, pos_end, i_start, i_end-1, f)
                    pos_start = pos_end
                else:
                    pos_start = offset
                i = i_end

        elif window_type == 'snp':
            # i: current line I'm reading
            # i_start: the position I aim to start a window
            # i_end: the position I aim to end a window

            i = 0
            i_start = offset

            while True:
                if i == n_snps: break

                chr = data.chromo[i]

                while True:
                    if i == i_start:
                        i_start = i
                        break
                    elif i == n_snps:
                        i_start = 'end'
                        break
                    elif data.chromo[i] != chr:
                        chr = data.chromo[i]
                        i_start = i + offset
                    else:
                        i += 1

                if i_start == 'end': break


                i = i_start
                chr = data.chromo[i_start]
                i_end = i_start + window_size

                while True:
                    if i == i_end:
                        i_end = i
                        status = True
                        break
                    elif i == n_snps:
                        i_end = 'end'
                        status = False
                        break
                    if data.chromo[i] != chr:
                        i_end = i
                        status = False
                        break
                    else:
                        i += 1

                if i_end == 'end': break
                if status:
                    calculate_statistics(
                        data.chromo[i_start], data.position[i_start], data.position[i_end-1],
                        i_start, i_end-1, f)
                    i_start = i_end
                else:
                    # i_start = offset  # this could have been the cause of the bug
                    i_start = i + offset
                i = i_end


os.system("cat {} | sort -V -k1,1 -k2,2 > {}".format(tmp_outfile, outfile))
os.system("rm {}".format(tmp_outfile))
