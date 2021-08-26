import numpy as np
from scipy.stats import entropy
from collections import Counter
import pandas as pd

def positional_entropy(sequences):
    entropy_vals = []
    for i in range(len(sequences[0])):
        nucs = [seq[i] for seq in sequences]
        ctr = Counter(nucs)
        vals =[v/len(nucs) for v in ctr.values()]
        entropy_vals.append(entropy(vals))
        
    return np.mean(entropy_vals)

def make_logomaker_df(sequences):
    df = pd.DataFrame()
    for i in range(len(sequences[0])):
        nucs = [seq[i] for seq in sequences]
        ctr = Counter(nucs)
        ctr ={k: v/len(nucs) for k, v in ctr.items()}
        df = df.append(ctr, ignore_index=True)

    df = df.fillna(0)
    return df
