import numpy as np
import sys


input_file = sys.argv[1]
output_file = sys.argv[2]
component = int(sys.argv[3])

h = 1000.0
#h = 1.0

in_mat = np.loadtxt(input_file, skiprows=1)
out_mat = in_mat
print(in_mat)
if component == 0:
    out_mat += h
elif component == 1:
    out_mat[0,0] += h
elif component == 2:
    out_mat[1,1] += h
elif component == 3:
    out_mat[2,2] += h
elif component == 12: # ancestral to 1 and 2
    out_mat[0,1] += h
    out_mat[1,0] += h
    out_mat[0,0] += h
    out_mat[1,1] += h
elif component == 23: # introgression from 2 to 3
    out_mat[1,2] += h
    out_mat[2,1] += h
print(out_mat)

np.savetxt(output_file, out_mat, header='{} {}'.format(in_mat.shape[0], in_mat.shape[1]), comments='')
