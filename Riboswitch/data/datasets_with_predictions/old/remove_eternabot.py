import pandas as pd

df1 = pd.read_json('subset_Z_est_Z_29Apr2020.json')
df2 = pd.read_json('subset_bp_est_bps_29Apr2020.json')

df1 = df1.loc[df1['Player'] != 'theeternabot']
df1.to_json('subset_est_Z_11May2020.json.zip')

df2 = df2.loc[df2['Player'] != 'theeternabot']
df2.to_json('subset_est_bps_11May2020.json.zip')

df3 = pd.read_json('all_packages_Z_est_Z_29Apr2020.json')
df4 = pd.read_json('all_packages_bp_est_bps_29Apr2020.json')

df3.to_json('all_packages_est_Z_11May2020.json.zip')
df4.to_json('all_packages_est_bps_11May2020.json.zip')
