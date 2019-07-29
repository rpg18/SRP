#Quick script to merge counts for duplicate genes produced by mmquant

import pandas as pd
df = pd.read_csv('All_Out_mmquant.csv', sep="\t", header = 0)
merge = df.groupby(['Gene']).agg('sum')	
merge.to_csv("Merged_out_mmquant.csv", sep='\t', header=True)
