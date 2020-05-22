import pandas as pd
import sys

df1 = pd.read_json(sys.argv[1])
df2 = pd.read_json(sys.argv[2])

cols_to_use = df2.columns.difference(df1.columns)
print(cols_to_use)
dfNew = pd.merge(df1, df2[cols_to_use], left_index=True, right_index=True, how='outer')
dfNew.to_json(sys.argv[3])
