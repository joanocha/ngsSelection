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

data = pd.read_csv(infile, '\t', names=["chromo", "position", "A01", "B01", "A02", "B02", "A12", "B12"])


a01 = data.A01
a02 = data.A02
a12 = data.A12

b01 = data.B01
b02 = data.B02
b12 = data.B12


n_snps = data.shape[0]

cum_a01 = np.cumsum(a01)
cum_a02 = np.cumsum(a02)
cum_a12 = np.cumsum(a12)
cum_b01 = np.cumsum(b01)
cum_b02 = np.cumsum(b02)
cum_b12 = np.cumsum(b12)




def calculate_statistics(chr, pos_start, pos_end, i_start, i_end, f):
    if i_start > i_end: return

    nsites= i_end - i_start + 1

    afd_01 = cum_a01[i_end] - (cum_a01[i_start-1] if i_start>0 else 0)
    afd_02 = cum_a02[i_end] - (cum_a02[i_start-1] if i_start>0 else 0)
    afd_12 = cum_a12[i_end] - (cum_a12[i_start-1] if i_start>0 else 0)

    total_afd_01 = cum_b01[i_end] - (cum_b01[i_start-1] if i_start>0 else 0)
    total_afd_02 = cum_b02[i_end] - (cum_b02[i_start-1] if i_start>0 else 0)
    total_afd_12 = cum_b12[i_end] - (cum_b12[i_start-1] if i_start>0 else 0)

    fst01_w = afd_01.sum() / total_afd_01.sum()
    fst02_w = afd_02.sum() / total_afd_02.sum()
    fst12_w = afd_12.sum() / total_afd_12.sum()

    fst01_uw = (afd_01 / total_afd_01).mean()
    fst02_uw = (afd_02 / total_afd_02).mean()
    fst12_uw = (afd_12 / total_afd_12).mean()

    def calculate_pbs(fst01, fst02, fst12):
        T1 = -np.log(1-fst01)
        T2 = -np.log(1-fst02)
        T3 = -np.log(1-fst12)
        pbs1 = (T1 + T2 - T3) / 2
        pbs2 = (T1 + T3 - T2) / 2
        pbs3 = (T2 + T3 - T1) / 2
        return pbs1, pbs2, pbs3

    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}\n'.format(
        chr, pos_start, pos_end,
        fst01_w, fst02_w, fst12_w,
        fst01_uw, fst02_uw, fst12_uw,
        *calculate_pbs(fst01_w, fst02_w, fst12_w), *calculate_pbs(fst01_uw, fst02_uw, fst12_uw),
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
