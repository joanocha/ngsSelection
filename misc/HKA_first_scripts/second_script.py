import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

window_size = 50000
step_size = 10000
offsets = np.arange(0, window_size, step_size)

threshold = 0.05
f_min = threshold
f_max = 1-threshold

data = pd.read_csv('vr_vv_inputHKA', '\t')
# data = pd.read_csv('input_small', '\t')

f_x = data.knownEM_x
f_y = data.knownEM_y

A = (f_x > f_min) & (f_x < f_max)
B = (f_y > f_min) & (f_x < f_max)
C = (f_x >= f_max) & (f_y <=f_min)
D = (f_y >= f_max) & (f_x <=f_min)

global_A = A.sum()
global_B = B.sum()
global_C = C.sum()
global_D = D.sum()


#global_A_AC_ratio = global_A / (global_A + global_C)
#global_B_BD_ratio = global_B / (global_B + global_D)

#print(global_A, global_B, global_C, global_D, global_A_AC_ratio, global_B_BD_ratio)
#raise Exception('debug')

n_snps = data.shape[0]
cum_A = np.cumsum(A)
cum_B = np.cumsum(B)
cum_C = np.cumsum(C)
cum_D = np.cumsum(D)


for offset in offsets:
    chr = data.chromo[0]
    chr_here = chr
    i = 0
    pos_start = offset
    pos_here = data.position[0]

    full_done = False
    while i < n_snps and not full_done:
        done = False
        while not done and not full_done:
            while pos_here < pos_start and chr_here == chr and i < n_snps-1 and not done:
                i += 1
                if i >= n_snps:
                    full_done = True
                    break
                pos_here = data.position[i]
                chr_here = data.chromo[i]
            if full_done: break
            if i >= n_snps-1:
                full_done = True
                break
            if chr_here == chr:
                i_start = i
                done = True
                break
            else:
                chr = chr_here
                pos_start = offset
        if full_done: break

        while chr_here == chr and pos_here-pos_start < window_size and i < n_snps-1:
            i += 1
            if i >= n_snps:
                full_done = True
                break
            chr_here = data.chromo[i]
            pos_here = data.position[i]
        if full_done: break

        i_end = i-1
        if i_end >= n_snps-1:
            full_done = True
            break

        if i_end >= i_start and chr == chr_here:
            a = cum_A[i_end] - (cum_A[i_start-1] if i_start>0 else 0)
            b = cum_B[i_end] - (cum_B[i_start-1] if i_start>0 else 0)
            c = cum_C[i_end] - (cum_C[i_start-1] if i_start>0 else 0)
            d = cum_D[i_end] - (cum_D[i_start-1] if i_start>0 else 0)

            if a + c > 0:
                HKA_x = -np.log10(chi2_contingency(np.array([[a, c], [global_A, global_C]]))[1])
            else:
                HKA_x = np.nan
            if b + d > 0:
                HKA_y = -np.log10(chi2_contingency(np.array([[b, d], [global_B, global_D]]))[1])
            else:
                HKA_y = np.nan
            if a + b > 0 and c + d > 0 and a + c > 0 and b + d > 0:
                HOMO = -np.log10(chi2_contingency(np.array([[a, c], [b, d]]))[1])
            else:
                HOMO = np.nan
            pos_end = pos_start + window_size
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chr, pos_start, pos_end, a, b, c, d, HKA_x, HKA_y, HOMO, c/a, d/b))

        if chr == chr_here:
            pos_start = pos_start + window_size
        else:
            pos_start = offset
        chr = chr_here
