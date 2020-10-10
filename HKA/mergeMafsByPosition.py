import pandas as pd
import numpy as np
import sys


inputs = sys.argv[1:-1]
output = sys.argv[-1]

print(inputs)

df = pd.read_csv(inputs[0], '\t', header=0,  usecols=["chromo", "position", "knownEM"])
df = df.rename(columns={'knownEM': "freq"+"_"+inputs[0].split('_')[0]})

print(df)

for fname in inputs[1:]:
    df2 = pd.read_csv(fname, '\t', header=0,  usecols=["chromo", "position", "knownEM"])
    df2 = df2.rename(columns={'knownEM': "freq"+"_"+fname.split('_')[0]})
    df = pd.merge(df, df2, how='inner', on=['chromo', 'position'])

print(df)
df.to_csv(output, sep='\t', index=False)
