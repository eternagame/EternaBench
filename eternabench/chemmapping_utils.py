import numpy as np
import pandas as pd
import sys
sys.path.append('/Users/hwayment/das/github/')

from RDATKit import rdatkit

import subprocess as sp
from arnie.bpps import bpps
from glob import glob
import sys, os
from copy import deepcopy

def write_rdat_as_dataframe(rdat_file, verbose=True):
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
            dct_to_add.update({'construct': c, 'filename': rdat.filename, 'seqpos': [x-1 for x in rdat.constructs[c].seqpos]})
            
            for k in annotations:
                if k in dct_to_add.keys():
                    if k=='modifier':
                        if not printed_modifier_warning:
                            print('Note: updating modifier for some')
                            printed_modifier_warning=True
                        dct_to_add.update({k: ', '.join(rdat.constructs[c].data[dat_ind].annotations[k][0])})
                        
                    elif k=='MAPseq':
                        # parse out sequence-specific MAPseq tags: like project name, design name, etc
                        for field in rdat.constructs[c].data[dat_ind].annotations[k]:
                            key = field.split(':')[0]
                            subfield = ':'.join(field.split(':')[1:])
                            dct_to_add[key] = subfield
                    else:
                        dct_to_add[k].append(rdat.constructs[c].data[dat_ind].annotations[k][0])
                        
                else:
                    dct_to_add.update({k: rdat.constructs[c].data[dat_ind].annotations[k][0]})

            dct_to_add.update({'reactivity': rdat.constructs[c].data[dat_ind].values,
                               'errors': rdat.constructs[c].data[dat_ind].errors})
            df = df.append(dct_to_add, ignore_index=True)
            
    return df

def filter_dataframe_with_cdhit(df, rdat_identifier):
    
    #write fasta file for CD-HIT-EST
    with open('%s.fasta' % rdat_identifier, 'w') as f:
        for i,seq in enumerate(df['sequence']):
            f.write(">%d\n" % i)
            f.write("%s\n" % seq[:-1*len('AAAAGAAACAACAACAACAAC')]) # removing tail
    print('Running CD-HIT-EST on %s...' % rdat_identifier)
    
    p = sp.Popen(('%s/cd-hit-est -i %s.fasta -o %s_cdhit_output -c 0.8' % (os.environ['CDHIT_PATH'], rdat_identifier,rdat_identifier)).split(' '),
                 stdout=sp.PIPE, stderr=sp.PIPE)
    p.wait()
    clusters, local_clust=[],[]    
    df['passed_CDHIT_filter'] = False
    
    print('Getting sequences from each cluster')
    with open('%s_cdhit_output.clstr' % rdat_identifier,'r') as f:
        for line in f.readlines():
            if not line.startswith('>'):
                if '>' in line:
                    local_clust.append(int(line.split('>')[1].split('.')[0]))
            else:
                clusters.append(local_clust)
                local_clust=[]

    for cluster in clusters:
        if len(cluster) > 0:
            if 'signal_to_noise' in df.keys():
                clust_stn = [float(x.split(':')[-1]) for x in df.iloc[cluster]['signal_to_noise'].values]
                df.loc[cluster[np.argmax(clust_stn)],['passed_CDHIT_filter']]=True
            else:
                df.loc[cluster[0],['passed_CDHIT_filter']]=True
                
    return df

def filter_FMN_containing(df):
    df['chem_string'] = [''.join(chem_list) for chem_list in df['chemical']]
    filtered_df = df.loc[~df['chem_string'].str.contains('FMN:[1-9]', regex=True)]
    filtered_df.drop(columns=['chem_string'],inplace=True)
    return filtered_df

def get_polyA_indicator(row, polyA_len=3, DEBUG=False):
    '''get indicator for which nucleotides are in polyA region.'''

    N = len(row['trimmed_sequence'])
    indicator=[0]*(polyA_len-1)

    for i in range(polyA_len,N+1):
        if ''.join(row['trimmed_sequence'][i-polyA_len:i]) == 'A'*polyA_len:
            indicator.append(1)
        else:
            indicator.append(0)

    if DEBUG:
        if len(indicator) != N:
            print(row['trimmed_sequence'])
            print(row['sequence'])
            print(indicator)

    return indicator
        
def write_concatenated_dataframe(df, reactivity_field='reactivity'):
    '''Input is dataframe with separate fields for each design, as well as possibly fields "p_<package_identifier>" 
    containing predicted p_unp vector for each design.
    Writes dataframe where all nucleotides are split up.
    '''
    
    # cut out parts of sequence that don't have reactivity data
    df['trimmed_sequence'] = df.apply(lambda row: [list(row['sequence'])[x] for x in row['seqpos'] if x < len(row['sequence'])], axis=1)
    n_constructs=len(df)
    
    df['length'] = [len(x) for x in df['sequence']]
    df = df.loc[df.length>6] # only use seqs longer than length that will be trimmed at polyA stage

    # doing this here to only apply to trimmed_sequence, but could be done earlier
    df['in_polyA'] = df.apply(lambda row: get_polyA_indicator(row, polyA_len=6), axis=1)

    # bad hack (and thats saying something for all the data handling hacks in here)
    # ... this is to get rid of reactivity values from constructs
    # where the length of the sequence is shorter than the reactivities recorded .. some of these in R75

    df[reactivity_field] = df.apply(lambda row: row[reactivity_field][:len(row['sequence'])], axis=1)


    concat_data = {reactivity_field: np.concatenate([x for x in df[reactivity_field]]),
                   'in_polyA': np.concatenate([x for x in df['in_polyA']]),
                  'nucleotide': np.concatenate([x for x in df['trimmed_sequence']])}
    
    if 'errors' in df.keys():
        df['errors'] = df.apply(lambda row: row['errors'][:len(row['sequence'])], axis=1)
        concat_data['errors'] = np.concatenate([x for x in df['errors']])
    
    #propagate construct-level information from before
    keys = [x for x in df.keys() if x not in [reactivity_field,'sequence','trimmed_sequence','seqpos','errors', 'in_polyA']]
    print('n_constructs',n_constructs)

    for key in keys:
        dat_to_add = []
        for i in df.index:
            dat_to_add.extend([df[key][i]]*len(df['trimmed_sequence'][i]))
            
        concat_data[key] = dat_to_add
    
    package_list = [x for x in df.keys() if x.startswith('p_')]
    
    for pkg in package_list:
        concat_data[pkg] = np.concatenate([x for x in df[pkg]])

    # for k, v in concat_data.items():
    #     print(k, len(v))

    return pd.DataFrame(data=concat_data)

def bootstrap_inds(len_item):
    return np.random.choice(range(len_item), len_item)
    
def filter_data(cdf, winsorize_cutoff=95):
    '''Filter concatenated dataframe to remove reactivity values below zero and above percentile cutoff.'''
    cdf_filt = pd.DataFrame()
    
    orig_len = len(cdf)
    
    cdf = cdf.loc[cdf['reactivity']>0]
    
    cutoff = np.percentile(cdf['reactivity'].values,winsorize_cutoff)
    cdf_filt = cdf.loc[cdf['reactivity']<cutoff]
        
    new_len = len(cdf_filt)
    
    print('%d of %d nucleotides (%.2f%%) removed, cutoff = %.2f' % (orig_len-new_len, orig_len, 100*(orig_len-new_len)/orig_len, cutoff))
    return cdf_filt
        


