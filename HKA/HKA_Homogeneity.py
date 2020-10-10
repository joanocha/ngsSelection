import numpy as np
import pandas as pd
import sys
import os
from scipy.stats import chi2_contingency, chisquare


infile = sys.argv[1]
outfile=sys.argv[2]
window_size = int(sys.argv[3])
step_size = int(sys.argv[4])
threshold=float(sys.argv[5])
window_type=sys.argv[6]
fixed_to = sys.argv[7]

if window_type not in ['snp', 'bp']:
    print("please specify if you want to do windows on snps or basepairs with option snp or bp")


tmp_outfile = outfile + '.unsorted'

offsets = np.arange(0, window_size, step_size)

np.set_printoptions(precision=4)

f_min = threshold
f_max = 1-threshold

data = pd.read_csv(infile, '\t')

f_x = data.knownEM_x
if "knownEM_y" not in data.columns:
    data["knownEM_y"] = np.nan
f_y = data.knownEM_y


A = ((f_x > f_min) & (f_x < f_max)).astype(int)
B = ((f_y > f_min) & (f_y < f_max)).astype(int)
if fixed_to == 'ingroup_and_outgroup':
    C = ((f_x >= f_max) & (f_y <=f_min)).astype(int)
elif fixed_to == 'outgroup':
    C= (f_x >= f_max).astype(int)
elif fixed_to == 'outgroup_or_outgroup':
    C=((f_x >= f_max) | ((f_y >= f_max) & (f_x<=f_min))).astype(int)
else:
    raise Exception("fixed to not implemented")

D = ((f_y >= f_max) & (f_x <=f_min)).astype(int)

global_A = A.sum()
global_B = B.sum()
global_C = C.sum()
global_D = D.sum()

global_prop_A = global_A / (global_A + global_C)
global_prop_C = 1 - global_prop_A

global_prop_B = global_B / (global_B + global_D)
global_prop_D = 1 - global_prop_B


print(global_A, global_B, global_C, global_D)


if fixed_to == 'ingroup_and_outgroup':
    total = A + B + C + D
elif fixed_to == 'outgroup':
    total = A + C
elif fixed_to == 'outgroup_or_outgroup':
    total = A + C

nsites_real = (total > 0).sum()
print(nsites_real, 'nsites_real')


n_snps = data.shape[0]
cum_A = np.cumsum(A)
cum_B = np.cumsum(B)
cum_C = np.cumsum(C)
cum_D = np.cumsum(D)


def calculate_statistics(chr, pos_start, pos_end, i_start, i_end, f):
    if i_start > i_end: return


    a = cum_A[i_end] - (cum_A[i_start-1] if i_start>0 else 0)
    b = cum_B[i_end] - (cum_B[i_start-1] if i_start>0 else 0)
    c = cum_C[i_end] - (cum_C[i_start-1] if i_start>0 else 0)
    d = cum_D[i_end] - (cum_D[i_start-1] if i_start>0 else 0)


    if a + c > 0:
        #HKA_x = -np.log10(chi2_contingency(np.array([[a, c], [global_A, global_C]]))[1]) #### old version
        n_sites = a + c
        expected_a = n_sites *  global_prop_A
        expected_c = n_sites * global_prop_C
        HKA_x_chi2, px  = chisquare([a, c], [expected_a, expected_c])
        HKA_x_logpval = -np.log10(px)
    else:
        HKA_x_chi2 = np.nan
        HKA_x_logpval = np.nan
    if b + d > 0:
        n_sites = b + d
        expected_b = n_sites *  global_prop_B
        expected_d = n_sites * global_prop_D
        HKA_y_chi2, py  = chisquare([b, d], [expected_b, expected_d])
        HKA_y_logpval = -np.log10(py)
        #HKA_y = -np.log10(chi2_contingency(np.array([[b, d], [global_B, global_D]]))[1])
        #chi2_y = chi2_contingency(np.array([[b, d], [global_B, global_D]]))[0]
    else:
        HKA_y_chi2 = np.nan
        HKA_y_logpval = np.nan
    if  a + b > 0 and c + d > 0 and a + c > 0 and b + d > 0:
        HOMO_chi2, ph, _, _  =  chi2_contingency(np.array([[a, c], [b, d]]))
        HOMO_logpval = -np.log10(ph)
    else:
        HOMO_chi2 = np.nan
        HOMO_logpval = np.nan

    f.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7:.4g}\t{8:.4g}\t{9:.4g}\t{10:.4g}\t{11:.4g}\t{12:.4g}\n'.format(chr, pos_start, pos_end, a, b, c, d, HKA_x_chi2, HKA_y_chi2, HOMO_chi2, HKA_x_logpval, HKA_y_logpval, HOMO_logpval))


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
                    i_start = i + offset
                i = i_end


os.system("cat {} | sort -V -k1,1 -k2,2 > {}".format(tmp_outfile, outfile))
os.system("rm {}".format(tmp_outfile))
