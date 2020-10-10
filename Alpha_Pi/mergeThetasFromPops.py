import pandas as pd
import numpy as np
import sys

pop1= sys.argv[1]
pop2= sys.argv[2]
pop3 = sys.argv[3]
output= sys.argv[4]

def load_data(pop):
   df=pd.read_csv(pop, '\t', header=None, comment='#', names=['information', 'chromo', 'wincenter',  'tP', 'nsites'])
   print(df.head())
   df['relative_pi'] = df['tP'] / df['nsites']
   df['winstart'] = df.information.str.split("(").str[3].str.split(")").str[0].str.split(",").str[0].astype(int)
   df['winend'] = df.information.str.split("(").str[3].str.split(")").str[0].str.split(",").str[1].astype(int)
   df = df.filter(items = ['chromo', 'winstart', 'winend', 'wincenter', 'tP', 'nsites', 'relative_pi'])
   print(df)
   return df

df1=load_data(pop1)
df2=load_data(pop2)
df3=load_data(pop3)

merged_pi = df1.merge(df2,  how='inner', on=["chromo", "winstart", "winend", "wincenter"])
merged_pi = merged_pi.merge(df3,  how='inner', on=["chromo", "winstart", "winend", "wincenter"])
merged_pi = merged_pi.rename(columns={'tP':'tP_z', 'nsites':'nsites_z', 'relative_pi':'relative_pi_z'})
print(merged_pi)
merged_pi['relative_total_pi'] = merged_pi['relative_pi_x'] + merged_pi['relative_pi_y'] + merged_pi['relative_pi_z']
merged_pi['relative_alpha'] = merged_pi['relative_pi_x'] / merged_pi['relative_total_pi']
merged_pi['total_pi'] = merged_pi['tP_x'] + merged_pi['tP_y'] + merged_pi['tP_z']
merged_pi['alpha'] = merged_pi['tP_x'] / merged_pi['total_pi']


print(merged_pi.head())

merged_pi.to_csv(output, sep='\t', index=False)
