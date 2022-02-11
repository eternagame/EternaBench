from polyleven import levenshtein
from tqdm import tqdm
import sys, os
import numpy as np


def window_seq(seq, window_size=100):
    n=len(seq)
    m=window_size
    frags=[]
    for i in range(0, n - m + 1):
        tmp = seq[i:i + m]
        frags.append(tmp)
        
    return frags
    
def empirically_determine_cutoff(seq_len):
    dists=[]
    for _ in range(100):
        seq1=''.join(np.random.choice(list('ACGU'), size=seq_len))
        seq2=''.join(np.random.choice(list('ACGU'), size=seq_len))
        dists.append(levenshtein(seq1, seq2)/seq_len)

    return np.percentile(dists, 0.025)

def compute_pairwise_with_windowing(seqs1, seqs2, MULT=1.1,length_minimum=70,cutoff=0.2):
    '''Return number of sequences in seqs2 that are within cutoff of any sequence in seqs1.'''
    count=0
    
    seqs1 = sorted(seqs1, key=len)
    
    within_range_seqs=[]

    for ind_i in tqdm(range(len(seqs2))):
        done=False
        for ind_j in range(len(seqs1)):
            
            shorter_seq, longer_seq = sorted([seqs1[ind_j], seqs2[ind_i]], key=len)

            if len(shorter_seq)>length_minimum and len(longer_seq) > MULT*len(shorter_seq):
                seq_frags = window_seq(longer_seq, window_size=len(shorter_seq)) # split into windows
                
            else:
                seq_frags = [longer_seq]

            for seq in seq_frags:
                dist = levenshtein(shorter_seq, seq)
                if dist/len(shorter_seq) < cutoff: #(1-cutoff):
                    # print(dist/len(shorter_seq))
                    # print('cutoff: ',cutoff_dct[len(shorter_seq)])
                    # print(shorter_seq)
                    # print(seq)
                    within_range_seqs.append(seqs2[ind_i])
                    done=True
                    break

            if done:
                #print('count', count)
                count+=1
                break
                    
    return count, within_range_seqs

def read_fasta(fasta):
    seqs=[]
    ids=[]
    with open(fasta, 'r') as f:
        for lin in f.readlines():
            if not lin.startswith('>'):
                seqs.append(lin.strip().upper())
            else:
                ids.append(lin.strip())
                
    return seqs, ids

def filter_fasta(fasta1, fasta2, outfile='OUT.fasta'):
    '''Return fasta of sequences in seqs2 that are nonredundant with the sequences in seqs1.'''

    print('Comparing %s (fil1) and %s (fil2)' % (fasta1, fasta2))
    seqs1, _ = read_fasta(fasta1)
    seqs2, ids2 = read_fasta(fasta2)

    count, within_range_seqs = compute_pairwise_with_windowing(seqs1, seqs2)
    
    print('%d sequences from %s are redundant with %s >50' % (count, fasta2, fasta1))
    print('Writing nonredundant seqs to %s' % outfile)
    
    with open(outfile,'w') as out:
        for seq in seqs2:
            if seq not in within_range_seqs:
                ind = ids2[seqs2.index(seq)]
                out.write(ind+'\n')
                out.write(seq+'\n')
                
    return count
                
fil1 = sys.argv[1]
fil2 = sys.argv[2]

outfile = os.path.basename(fil1)+'_'+os.path.basename(fil2)
outfile = '/home/groups/rhiju/hannahw1/EternaBench/CrossDatasetAnalysis/lev_output/'+outfile.replace('.fasta','')
count = filter_fasta(fil1, fil2, outfile=outfile+'80_windowed.fasta')
            
