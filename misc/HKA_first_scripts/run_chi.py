import numpy as np
import pandas as pd
from scipy.stats import chisquare


data = pd.read_csv('result.txt', '\t', header=None, names=['chr', 'start', 'end', 'f', 'p', 'ratio'])
global_ratio = data.f.sum() / (data.f.sum()+data.p.sum())

n_windows = data.shape[0]

p_values = np.empty(n_windows)

print(n_windows)

for i in range(n_windows):
    if i % 10000 == 0: print(i)
    p = data.p[i]
    f = data.f[i]
    expected = np.array([global_ratio, 1-global_ratio]) * (f+p)
    observed = [f, p]
    p_values[i] = chisquare(observed, expected)[1]


np.savetxt('p_values.txt', p_values)
