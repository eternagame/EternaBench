import numpy as np
import pandas as pd
import subprocess as sp
import arnie.utils as utils

def get_lig_constraint(row):
    seq = row['Sequence']
    lig = row['ligand']
    
    try:

        lig_aptamer = utils.write_constraints(seq, LIG=True, 
            lig1 = aptamers[lig][0],lig2 = aptamers[lig][1])

        if len(lig_aptamer) > len(seq):
                lig_aptamer = utils.write_constraints(seq, LIG=True,
                    lig1 = aptamers['%s_rev' % lig][0],lig2 = aptamers['%s_rev' % lig][1])
    except:
        lig_aptamer='.'
        print(lig_aptamer)


    #print(seq)

    return lig_aptamer

def get_MS2_constraint(row):
    try:
        constr = utils.write_constraints(row['Sequence'], MS2=True)
    except:
        constr = '.'
        print(constr)
    return constr
        

def get_lig_MS2_constraint(row):
    try:
        seq= row['Sequence']
        lig = row['ligand']
        lig_ms2_aptamer = utils.write_constraints(seq, MS2=True, LIG=True, lig1 = aptamers[lig][0],lig2 = aptamers[lig][1])

        if len(lig_ms2_aptamer) > len(seq):
            lig_ms2_aptamer = utils.write_constraints(seq, MS2=True, LIG=True, lig1 = aptamers['%s_rev' % lig][0],lig2 = aptamers['%s_rev' % lig][1])
    except:
        lig_ms2_aptamer='.'

    return lig_ms2_aptamer

fields_w_cluster = ['Design', 'Sequence', 'Player', 'Puzzle_Name', 'FoldChange','KDON', 'KDOFF', 'NumberOfClusters']
fields = ['Design','Sequence', 'Player', 'Puzzle_Name','FoldChange', 'KDON', 'KDOFF', 'Puzzle_Type']

tmp_df_list=[]
for rnd in [93,95,96,97,98]:
    tmp_df = pd.read_csv('R%d_fmn.txt' % rnd, delimiter='\t')[fields_w_cluster]
    tmp_df = tmp_df.loc[tmp_df['NumberOfClusters']>50]
    tmp_df['Round'] = rnd
    tmp_df['ligand'] = 'FMN'
    tmp_df['min_kd_val'] = np.min([np.min(tmp_df['KDON']),np.min(tmp_df['KDOFF'])])
    tmp_df_list.append(tmp_df)
    
df_101 = pd.read_csv('R101.txt', delimiter='\t')[fields]
df_101 = df_101.loc[df_101['Puzzle_Type']=='FMN']
df_101['ligand'] = 'FMN'
df_101['Round'] = 101
df_101 = df_101.drop(columns=['Puzzle_Type'])
df_101['min_kd_val'] = np.min([np.min(df_101['KDON']),np.min(df_101['KDOFF'])])

tmp_df_list.extend([df_101])

df_107 = pd.read_csv('R107.txt',delimiter='\t')[fields]
df_107 = df_107[df_107.Puzzle_Type.notnull()]
df_107['ligand'] = df_107.apply(lambda row: row['Puzzle_Type'].split('_')[0],axis=1)
df_107['output'] = df_107.apply(lambda row: row['Puzzle_Type'].split('_')[1],axis=1)
df_107 = df_107.loc[df_107['output']=='MS2']
df_107 = df_107.drop(columns=['output','Puzzle_Type'])
df_107['min_kd_val'] = np.min([np.min(df_107['KDON']),np.min(df_107['KDOFF'])])
df_107['Round'] = 107

tmp_df_list.extend([df_107])


df = pd.concat(tmp_df_list,ignore_index=True)
 
df.loc[df['Puzzle_Name'].str.contains('Same State'),'switch'] = 'ON'
df.loc[df['Puzzle_Name'].str.contains('Exclusion'),'switch']= 'OFF'

df = df.rename(columns={'FoldChange': 'Activation Ratio'})
df = df.dropna(subset=['KDON','KDOFF'])

#filter nonstandard motifs
aptamers = {'FMN':[('nAGGAUAU', '(xxxxxx('),('AGAAGGn', ')xxxxx)')],
            'FMN_rev':[('AGAAGGn', '(xxxxx('),('nAGGAUAU', ')xxxxxx)')],
            'Theo':[('GAUACCAG','(xxx(((('),('CCCUUGGCAGC',')xxx)))xxx)')],
            'Theo_rev':[('CCCUUGGCAGC','(xxx(((xxx('),('GAUACCAG',')xxx))))')],
            'Trp':[('AGGACCGG','((xxx((('),('CCGCCACU',')))xxx))')],
            'Trp_rev':[('CCGCCACU','(((xxx(('),('AGGACCGG','))xxx)))')]}

#MS2:
df = df.loc[df['Sequence'].str.contains('ACAUGAGGAUCACCCAUGU')]

tmp_df_list=[]
for lig in ['FMN','Trp','Theo']:
    tmp_df = df.loc[df['ligand']==lig]
    tmp_df = tmp_df.loc[tmp_df['Sequence'].str.contains(aptamers[lig][0][0].replace('n',''))]
    tmp_df_list.extend([tmp_df])

df = pd.concat(tmp_df_list,ignore_index=True)

df['MS2_aptamer'] = df.apply( lambda row: get_MS2_constraint(row), axis=1)
df['lig_aptamer'] = df.apply( lambda row: get_lig_constraint(row), axis=1)
df['MS2_lig_aptamer'] = df.apply( lambda row: get_lig_MS2_constraint(row),axis=1)
df['constraints_worked'] = df.apply( lambda row: len(row['Sequence']) - len(row['MS2_lig_aptamer']), axis=1)
df = df.loc[df['constraints_worked']==0]
df = df.dropna(subset=['Sequence'])

df['logkd_nolig'] = np.log(np.where(df['switch']=='OFF', df['KDON'], df['KDOFF']))
df['logkd_lig'] = np.log(np.where(df['switch']=='OFF', df['KDOFF'], df['KDON']))
df['logkd_nolig_scaled'] = df['logkd_nolig'] - np.log(df['min_kd_val']) # log min value in dataset, k^*
df['logkd_lig_scaled'] = df['logkd_lig'] - np.log(df['min_kd_val']) # log min value in dataset, k^*

with open('eterna_switches.fasta', 'w') as f:
    for i,seq in enumerate(df['Sequence'].values):
        f.write(">%d\n" % i)
        f.write("%s\n" % seq)
print('Running CD-HIT-EST')

p = sp.Popen(('/Users/hwayment/das/github/cdhit/cd-hit-est -i eterna_switches.fasta -o eterna_cdhit_output -c 0.8').split(' '),
             stdout=sp.PIPE, stderr=sp.PIPE)
clusters=[]
local_clust=[]

p.wait()

with open('eterna_cdhit_output.clstr','r') as f:
    for line in f.readlines():
        if not line.startswith('>'):
            if '>' in line:
                local_clust.append(int(line.split('>')[1].split('.')[0]))
        else:
            clusters.append(local_clust)
            local_clust=[]
            
filtered_inds=[]
for x in clusters:
    if len(x)>0:
        filtered_inds.append(x[0])

df = df.loc[filtered_inds]
print(len(df))
df['sequence'] = df['Sequence']

df.to_json('processed_eterna_df.json')
