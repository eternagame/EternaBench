import numpy as np
import pandas as pd
import subprocess as sp
import arnie.utils as utils
from datetime import date

todaysdate = date.today().strftime("%d%b%Y")

from pandarallel import pandarallel
pandarallel.initialize()

Round_dict={88: 1, 93:2, 95:3, 96:4, 97:5, 98:6, 101: 7, 107:8,'Ribologic':'Ribologic'}

def get_lig_constraint(row):
    seq = row['sequence']
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
        constr = utils.write_constraints(row['sequence'], MS2=True)
    except:
        constr = '.'
        print(constr)
    return constr 

def get_lig_MS2_constraint(row):
    try:
        seq= row['sequence']
        lig = row['ligand']
        lig_ms2_aptamer = utils.write_constraints(seq, MS2=True, LIG=True, lig1 = aptamers[lig][0],lig2 = aptamers[lig][1])

        if len(lig_ms2_aptamer) > len(seq):
            lig_ms2_aptamer = utils.write_constraints(seq, MS2=True, LIG=True, lig1 = aptamers['%s_rev' % lig][0],lig2 = aptamers['%s_rev' % lig][1])
    except:
        lig_ms2_aptamer='.'

    return lig_ms2_aptamer

fields_w_cluster = ['Design', 'Sequence', 'Player', 'Puzzle_Name', 'FoldChange','KDON', 'KDOFF', 'NumberOfClusters']
fields = ['Design','Sequence', 'Player', 'Puzzle_Name','FoldChange', 'KDON', 'KDOFF', 'Puzzle_Type']
r88_fields=['Design','Sequence', 'Player', 'Puzzle_Name','FoldChange', 'KDON', 'KDOFF',]
tmp_df_list=[]
for rnd in [88, 93,95,96,97,98]:
    print('#####', rnd)
    try:
        tmp_df = pd.read_csv(os.environ['ETERNABENCH_PATH']+'/data/RiboswitchPreprocessing/R%d_fmn.txt' % rnd, delimiter='\t')[fields_w_cluster]
        tmp_df = tmp_df.loc[tmp_df['NumberOfClusters']>50]

    except:
        tmp_df = pd.read_csv('R%d_fmn.txt' % rnd, delimiter='\t')[r88_fields]

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
df = df.rename(columns={'Sequence':'sequence'})

#####
### Add ribologic dataset

rdf = pd.read_csv(os.environ['ETERNABENCH_PATH']+'/data/RiboswitchPreprocessing/ribologic_SI.txt',delimiter='\t')

rdf = rdf.drop(['secstruct_ligand','secstruct_noligand'],axis=1)
rdf = rdf.loc[rdf["ligand"] != 'miRNA']
print(rdf.keys())
rdf['KDON'] = rdf['Kd_ON']
rdf['KDOFF'] = rdf['Kd_OFF']
rdf = rdf.dropna(subset=['sequence','Kd_ON','Kd_OFF'])
# rdf['logkd_nolig'] = np.log(np.where(rdf['switch']=='OFF', rdf['Kd_ON'], rdf['Kd_OFF']))
# rdf['logkd_lig'] = np.log(np.where(rdf['switch']=='OFF', rdf['Kd_OFF'], rdf['Kd_ON']))
# rdf['logkd_nolig_scaled'] = rdf['logkd_nolig'] - np.log(2) # log min value in dataset, k^*
# rdf['logkd_lig_scaled'] = rdf['logkd_lig'] - np.log(2) # log min value in dataset, k^*
# rdf['min_kd_val'] = 2

rdf['min_kd_val'] = np.min([np.min(rdf['KDON']),np.min(rdf['KDOFF'])])
rdf['Round'] = 'Ribologic'
rdf = rdf.rename(columns={'activation_ratio':'Activation Ratio'})
rdf.replace('Theophylline','Theo',inplace=True)
rdf.replace('Tryptophan','Trp',inplace=True)

rdf = rdf.reset_index()
df = df.append(rdf, ignore_index=True)
df = df.rename(columns={'Sequence':'sequence'})
df = df.dropna(subset=['sequence'])

#####


#filter nonstandard motifs
aptamers = {'FMN':[('nAGGAUAU', '(xxxxxx('),('AGAAGGn', ')xxxxx)')],
            'FMN_rev':[('AGAAGGn', '(xxxxx('),('nAGGAUAU', ')xxxxxx)')],
            'Theo':[('GAUACCAG','(xxx(((('),('CCCUUGGCAGC',')xxx)))xxx)')],
            'Theo_rev':[('CCCUUGGCAGC','(xxx(((xxx('),('GAUACCAG',')xxx))))')],
            'Trp':[('AGGACCGG','((xxx((('),('CCGCCACU',')))xxx))')],
            'Trp_rev':[('CCGCCACU','(((xxx(('),('AGGACCGG','))xxx)))')]}

#MS2:
print('Prior filtering MS2, ', len(df))
df = df.loc[df['sequence'].str.contains('ACAUGAGGAUCACCCAUGU')]
print('Post filtering MS2, ', len(df))

tmp_df_list=[]
for lig in ['FMN','Trp','Theo']:
    tmp_df = df.loc[df['ligand']==lig]
    tmp_df = tmp_df.loc[tmp_df['sequence'].str.contains(aptamers[lig][0][0].replace('n',''))]
    tmp_df_list.extend([tmp_df])

df = pd.concat(tmp_df_list,ignore_index=True)

print('writing constraints')
df['MS2_aptamer'] = df.parallel_apply( lambda row: get_MS2_constraint(row), axis=1)
df['lig_aptamer'] = df.parallel_apply( lambda row: get_lig_constraint(row), axis=1)
df['MS2_lig_aptamer'] = df.parallel_apply( lambda row: get_lig_MS2_constraint(row),axis=1)
df['constraints_worked'] = df.parallel_apply( lambda row: len(row['sequence']) - len(row['MS2_lig_aptamer']), axis=1)

print('wrote constraints')
df = df.loc[df['constraints_worked']==0]

df['logkd_nolig'] = np.log(np.where(df['switch']=='OFF', df['KDON'], df['KDOFF']))
df['logkd_lig'] = np.log(np.where(df['switch']=='OFF', df['KDOFF'], df['KDON']))

df['logkd_nolig_scaled'] = df['logkd_nolig'] - np.log(df['min_kd_val']) # log min value in dataset, k^*
df['logkd_lig_scaled'] = df['logkd_lig'] - np.log(df['min_kd_val']) # log min value in dataset, k^*

df['Dataset'] = df.apply(lambda row: str(Round_dict[row['Round']])+'_'+row['ligand'], axis=1 )

df = df.sort_values('Dataset')

print(df.groupby('Dataset').size())

# Reorder to put ribologic first
tmp_df1 = df.loc[df.Dataset.str.contains('Ribologic')]

tmp_df2 = df.loc[~df.Dataset.str.contains('Ribologic')]

df = tmp_df1.append(tmp_df2, ignore_index=True)

df = df.drop_duplicates(subset='sequence', keep='first')
print(df.groupby('Dataset').size())

df.to_json('unfiltered_dataset_06Aug2021.json.zip')

new_df = pd.DataFrame()

for rnd in df.Dataset.unique():
    tmp = df.loc[df.Dataset==rnd]
    tmp = tmp.reset_index()
    print(len(tmp))
    with open('eterna_switches_%s.fasta' % rnd, 'w') as f:
        for i,seq in enumerate(tmp['sequence'].values):
            f.write(">%d\n" % i)
            f.write("%s\n" % seq)
    print('Running CD-HIT-EST', rnd,len(tmp))

    p = sp.Popen(('%s/cdhit/cd-hit-est -i eterna_switches_%s.fasta -o eterna_cdhit_output_%s -c 0.8' % (os.environ['CDHIT_PATH'], rnd,rnd)).split(' '),
                 stdout=sp.PIPE, stderr=sp.PIPE)
    clusters=[]
    local_clust=[]

    p.wait()

    with open('eterna_cdhit_output_%s.clstr' % rnd,'r') as f:
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

    tmp['passed_CDHIT_filter'] = False
    tmp.loc[filtered_inds,'passed_CDHIT_filter'] = True
    new_df = new_df.append(tmp, ignore_index=True)

print(new_df.groupby(['Dataset', 'passed_CDHIT_filter']).size())

new_df['log_AR'] = np.log(new_df['Activation Ratio'])

new_df = new_df.drop(columns=['Kd_OFF_predicted','Kd_ON_predicted','AR_predicted','level_0'])
new_df.to_json(os.environ['ETERNABENCH_PATH']+'/data/EternaBench_Riboswitch_Full_%s.json.zip' % todaysdate)
new_df.loc[new_df.passed_CDHIT_filter].to_json(os.environ['ETERNABENCH_PATH']+'/data/EternaBench_Riboswitch_Filtered_%s.json.zip' % todaysdate)
