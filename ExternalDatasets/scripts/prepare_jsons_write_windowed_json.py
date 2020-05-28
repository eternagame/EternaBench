import numpy as np
import pandas as pd

def write_windowed_df(old_df, window_size=600, overlap=25):
    df = pd.DataFrame()
    for row in old_df.iterrows():
        seq_len = len(row[1]['sequence'])
        start_inds = [x for x in range(0, seq_len, window_size - overlap)]
        end_inds = [np.min([seq_len, x+window_size]) for x in start_inds]
        for s,e in list(zip(start_inds, end_inds)):
            seq_frag = row[1]['sequence'][s:e]
            reac_frag = row[1]['reactivity'][s:e]
            seqpos = [x for x in range(len(seq_frag))]
            df = df.append({'sequence': seq_frag, 'reactivity': reac_frag, 'seqpos': seqpos,
                            'ID': row[1]['ID'], 'info': row[1]['info']},ignore_index=True)
    return df

mod_df = pd.read_json('inputs/EB_parsed_nowindow_moderna_data.json')
mod_df['ID']=mod_df['name']
mod_df['info'] = ['Mauger,2019 ' + x + ' nLUC mRNA' for x in mod_df['name']]
mod_df.drop(columns='name',inplace=True)
polyA_df = pd.read_json('inputs/other_polyA_data_partial_df.json')
simon_df = pd.read_json('inputs/simon_partial_df.json')
polyA_df['sequence'] = polyA_df['seq']
polyA_df['reactivity'] = polyA_df['data']
df = pd.concat([mod_df, polyA_df, simon_df],ignore_index=True)
df = df.loc[~df['ID'].str.contains('denatured')]
df = df.loc[~df['ID'].str.contains('rep2')]
df.to_json('full_df.json')
df['seqpos'] = [[x for x in range(len(seq))] for seq in df['sequence']]
df['length'] = [len(x) for x in df['sequence']]

for w in [300,600,900]:
	wdf = write_windowed_df(df,window_size=w)
	wdf.to_json('full_df_window%d.json' % w)
