import pandas as pd
import sys
import numpy as np
infile=sys.argv[1]
outfile=sys.argv[2]

# unlike rule rank_genes this script considers the windows that have mover than one gene in the same window

df = pd.read_csv(infile, '\t')
print(df)
df = df[df["overlapping genes"] != "."]
df["overlapping genes"] =  df["overlapping genes"].str.split(",")
genes = pd.unique(np.concatenate(df["overlapping genes"].values))
pd.DataFrame(genes).to_csv(outfile, index=False, header=False)
