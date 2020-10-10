import numpy as np
import pandas as pd

threshold = 0.05
F_min = threshold
F_max = 1-threshold

data = pd.read_csv('test.txt', '\t', index_col=0)

Fx = data.knownEM_x
Fy = data.knownEM_y
Mx = data.major_x
mx = data.minor_x
My = data.major_y
my = data.minor_y

window_size = 50000

data['polymorphic'] = (Fx > F_min) & (Fx < F_max)

data['fixed_difference'] = (
    ((Fx<=F_min) & (Fy<=F_min) & (Mx!=My))
    |
    ((Fx<=F_min) & (Fy>=F_max) & (Mx!=my))
    |
    ((Fx>=F_max) & (Fy<=F_min) & (mx!=My))
    |
    ((Fx>=F_max) & (Fy>=F_max) & (mx!=my))
)

global_fixed = data.fixed_difference.sum()
global_polymorphic = data.polymorphic.sum()
global_ratio = global_fixed / (global_fixed+global_polymorphic)

n_snps = data.shape[0]
cumsum_p = data.polymorphic.cumsum().values
cumsum_f = data.fixed_difference.cumsum().values

chr = data.chromo[0]
chr_here = chr
i_start = 0
i = 0
pos_start = data.position[0]
pos_here = pos_start

while i < n_snps:
    while chr_here == chr and pos_here-pos_start < window_size:
        chr_here = data.chromo[i]
        pos_here = data.position[i]
        i += 1
    i_end = i-1
    f = cumsum_f[i_end] - cumsum_f[i_start]
    p = cumsum_p[i_end] - cumsum_p[i_start]
    ratio = f / (f + p)
    pos_end = data.position[i_end]
    print('{}\t{}\t{}\t{}\t{}\t{}'.format(chr, pos_start, pos_end, f, p, ratio))
    chr = chr_here
    i_start = i
    pos_start = pos_here
