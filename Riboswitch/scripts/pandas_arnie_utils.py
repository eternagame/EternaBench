import numpy as np
import pandas as pd
import os
from arnie.bpps import bpps
from arnie.pfunc import pfunc

def compute_bp_kfold_ref(**kwargs):
	return bpps('ACAUGAGGAUCACCCAUGU', **kwargs)[0][-1]

def compute_bp_kfold_ref_with_flanking(**kwargs):
    return bpps('GGGUAUGUCGCAGAAACAUGAGGAUCACCCAUGUAACUGCGACAUACCC', **kwargs)[15][-16]

def compute_Z_kfold_ref(**kwargs):

    Z_MS2 = pfunc('ACAUGAGGAUCACCCAUGU', constraint='(((((x((xxxx)))))))', **kwargs)
    Z = pfunc('ACAUGAGGAUCACCCAUGU', **kwargs)
    return np.log(Z/Z_MS2)

def compute_Z_kfold_ref_with_flanking(**kwargs):
    Z_MS2 = pfunc('GGGUAUGUCGCAGAAACAUGAGGAUCACCCAUGUAACUGCGACAUACCC',
     constraint='...............(((((x((xxxx)))))))...............', **kwargs)
    Z = pfunc('GGGUAUGUCGCAGAAACAUGAGGAUCACCCAUGUAACUGCGACAUACCC', **kwargs)
    return np.log(Z/Z_MS2)

def compute_bp_kfold(row, flanking=False, **kwargs):
    if flanking:
        return compute_bp_kfold_with_flanking(row, **kwargs)
        
    bp_matrix = bpps(row['sequence'], **kwargs)
    first_base = row['MS2_aptamer'].find('(')
    second_base = row['MS2_aptamer'].rfind(')')
    return bp_matrix[first_base][second_base]

def compute_bp_kfold_with_flanking(row, **kwargs):

    flank1='GGGUAUGUCGCAGAA'
    flank2='AACUGCGACAUACCC'
    seq = flank1+row['sequence']+flank2
    mod_ms2_aptamer = '.'*15+row['MS2_aptamer']+'.'*15
    bp_matrix = bpps(seq, **kwargs)
    first_base = mod_ms2_aptamer.find('(')
    second_base = mod_ms2_aptamer.rfind(')')
    return bp_matrix[first_base][second_base]

def compute_Z_kfold_with_flanking(row, **kwargs):

    flank1='GGGUAUGUCGCAGAA'
    flank2='AACUGCGACAUACCC'
    seq = flank1+row['sequence']+flank2
    mod_ms2_aptamer = '.'*15+row['MS2_aptamer']+'.'*15
    mod_lig_aptamer = '.'*15+row['lig_aptamer']+'.'*15
    mod_ms2_lig_aptamer = '.'*15+row['MS2_lig_aptamer']+'.'*15

    concentration = {'FMN':200e-6,'Theo':2e-3,'Trp':2.4e-3,'Theophylline':2e-3,'Tryptophan':2.4e-3}
    intrinsic_kd = {'FMN':2.2e-6,'Theo':20e-6,'Trp':13e-6, 'Theophylline': 20e-6, 'Tryptophan':13e-6}

    ligand_bonus = concentration[row['ligand']]/intrinsic_kd[row['ligand']]

    Z = pfunc(seq, **kwargs)
    Z_MS2 = pfunc(seq, constraint=mod_ms2_aptamer, **kwargs)
    Z_lig = pfunc(seq, constraint=mod_lig_aptamer, **kwargs)
    try:
        Z_MS2_lig = pfunc(seq, constraint=mod_ms2_lig_aptamer, **kwargs)
    except:
        Z_MS2_lig = 0

    return np.log((Z+ligand_bonus*Z_lig)/(Z_MS2+ligand_bonus*Z_MS2_lig)), np.log(Z/Z_MS2)

def compute_Z_kfold(row, flanking=False, **kwargs):
    if flanking:
        return compute_Z_kfold_with_flanking(row, **kwargs)
        
    concentration = {'FMN':200e-6,'Theo':2e-3,'Trp':2.4e-3,'Theophylline':2e-3,'Tryptophan':2.4e-3}
    intrinsic_kd = {'FMN':2.2e-6,'Theo':20e-6,'Trp':13e-6, 'Theophylline': 20e-6, 'Tryptophan':13e-6}

    seq=row['sequence']
    ligand_bonus = concentration[row['ligand']]/intrinsic_kd[row['ligand']]

    Z = pfunc(seq, **kwargs)
    Z_MS2 = pfunc(seq, constraint=row['MS2_aptamer'], **kwargs)
    Z_lig = pfunc(seq, constraint=row['lig_aptamer'], **kwargs)
    try:
        Z_MS2_lig = pfunc(seq, constraint=row['MS2_lig_aptamer'], **kwargs)
    except:
        Z_MS2_lig = 0

    return np.log((Z+ligand_bonus*Z_lig)/(Z_MS2+ligand_bonus*Z_MS2_lig)), np.log(Z/Z_MS2)

