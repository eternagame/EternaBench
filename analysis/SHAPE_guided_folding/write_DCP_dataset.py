import numpy as np
import pandas as pd

from glob import glob
import os, sys, pickle, requests

import chemmapping_utils as utils
from RDATKit import rdatkit
from copy import deepcopy

def write_rdat_as_dataframe_loc(rdat_file):
    df = pd.DataFrame()
    rdat = rdatkit.RDATFile()
    rdat.load(open(rdat_file))
    print(rdat.comments)
    printed_modifier_warning = False
    print(rdat.annotations)
    constructs = [x for x in rdat.constructs.keys()]
    for c in constructs:
        for dat_ind in range(len(rdat.constructs[c].data)):
            annotations = [x for x in rdat.constructs[c].data[dat_ind].annotations.keys()]
            dct_to_add = deepcopy(rdat.annotations)
            dct_to_add.update({'construct': c, 'filename': rdat.filename,
                               'sequence': rdat.constructs[c].sequence,
                               'seqpos': [x-1 for x in rdat.constructs[c].seqpos]})
            
            for k in annotations:
                if k in dct_to_add.keys():
                    if k=='modifier':
                        if not printed_modifier_warning:
                            print('Note: updating modifier for some')
                            printed_modifier_warning=True
                        dct_to_add.update({k: rdat.constructs[c].data[dat_ind].annotations[k][0]})
                        
                    elif k=='MAPseq':
                        # parse out sequence-specific MAPseq tags: like project name, design name, etc
                        for field in rdat.constructs[c].data[dat_ind].annotations[k]:
                            key, subfield = field.split(':')
                            dct_to_add[key] = subfield
                    else:
                        dct_to_add[k].append(rdat.constructs[c].data[dat_ind].annotations[k][0])
                        
                else:
                    dct_to_add.update({k: rdat.constructs[c].data[dat_ind].annotations[k][0]})
            dct_to_add.update({'reactivity': rdat.constructs[c].data[dat_ind].values,
                               'errors': rdat.constructs[c].data[dat_ind].errors})
            df = df.append(dct_to_add, ignore_index=True)
            
    return df
    
all_rdats = glob('raw_RDATS/DCP_RDATS_SHAPE_DMS/*.rdat')
print(all_rdats)

full_df = pd.DataFrame()
for rd in all_rdats:
    df = write_rdat_as_dataframe_loc(rd)
    full_df = full_df.append(df,ignore_index=True)
    
tmp = full_df.loc[full_df.modifier=='DMS']
tmp2 = pd.DataFrame()
for _, row in tmp.iterrows():
    if len(row['chemical'])==2:
        tmp2 = tmp2.append(row, ignore_index=True)

print('length', len(tmp2))
tmp2.to_json('ex_dataset.json')

