import pandas as pd
import os

df = pd.read_json('/Users/hwayment/das/github/EternaBench/data/EternaBench_ChemMapping_Filtered_10Jul2021.json.zip')
tmp = df.loc[df.Dataset=='Round 01']
tmp.to_json('/Users/hwayment/das/github/EternaBench/data/DEMO_ChemMapping.json.zip')


df = pd.read_json('/Users/hwayment/das/github/EternaBench/data/EternaBench_Riboswitch_Filtered_07Aug2021.json.zip')
tmp = df.loc[df.Dataset=='Ribologic_FMN']
tmp.to_json('/Users/hwayment/das/github/EternaBench/data/DEMO_Riboswitch.json.zip')

