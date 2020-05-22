import numpy as np
import pandas as pd
import os, argparse
from arnie.bpps import bpps
from pandarallel import pandarallel
import pandas_arnie_utils as utils
from datetime import date

todaysdate = date.today().strftime("%d%b%Y")

ARNIEDIR = '/home/users/hannahw1/arnie'

DIR = '/home/users/hannahw1/efold_params'

packages_all_options = {
'vienna_2': {'package': 'vienna_2'},
'vienna_2_nodangles': {'package':'vienna_2', 'dangles': False},
'vienna_2_60C': {'package': 'vienna_2', 'T': 60},
'vienna_1': {'package': 'vienna_1'},

 'nupack_99': {'package':'nupack'},
 'nupack_99_nodangles':{'package':'nupack', 'dangles': False},
 'nupack_95':{'package':'nupack_95'},
 'nupack_95_nodangles':{'package':'nupack_95', 'dangles': False},

'rnastructure': {'package':'rnastructure'},
'rnastructure_nocoax': {'package':'rnastructure', 'coaxial':False},

'contrafold_1': {'package':'contrafold_1'},
'contrafold_2': {'package':'contrafold_2'},
'contrafold_2_nc': {'package':'contrafold_2','param_file':'%s/parameter_files/contrafold.params.noncomplementary' % ARNIEDIR},
'learntofold': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/learntofold.contrafold.params' % ARNIEDIR},
'cyclefold': {'package': 'cyclefold'},

'vienna_langdon_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_langdon2018.par' % ARNIEDIR},
'vienna_rnasoft_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_andronescu2007.par' % ARNIEDIR},
'rnasoft_99' : {'package':'rnasoft_99'},
'rnasoft_07' : {'package':'rnasoft_07'},
'rnasoft_blstar' : {'package':'rnasoft'},
'rnasoft_99_nodangles' : {'package':'rnasoft_99-no-dangles'},
'rnasoft_bl_nodangles' : {'package':'rnasoft_bl-no-dangles'},
'rnasoft_lam-cg' : {'package':'rnasoft_lam-cg'},
'rnasoft_nom-cg' : {'package':'rnasoft_nom-cg'},

'eternafold_A': {'package': 'contrafold_2', 'param_file':'%s/params_A' % DIR},
'eternafold_B': {'package': 'contrafold_2', 'param_file':'%s/params_B' % DIR},
'eternafold_C': {'package': 'contrafold_2', 'param_file':'%s/params_C' % DIR},
'eternafold_D': {'package': 'contrafold_2', 'param_file':'%s/params_D' % DIR},
'eternafold_E': {'package': 'contrafold_2', 'param_file':'%s/params_E' % DIR},
'eternafold_F': {'package': 'contrafold_2', 'param_file':'%s/params_F' % DIR},
'eternafold_G': {'package': 'contrafold_2', 'param_file':'%s/params_G' % DIR},
 
 }

subset_list = ['vienna_2', 'nupack_99', 'rnastructure', 'contrafold_2', 'rnasoft_blstar', 'eternafold_B']

packages_Z_subset = {k: packages_all_options[k] for k in packages_all_options.keys()
     if k.startswith('vienna') or k.startswith('contrafold') or k.startswith('eternafold') or k.startswith('rnastructure')}

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Estimate protein affinity for sequences in input dataset.
    Input: .json file (eterna or ribologic dataset)
    Output: .json file (same as input, but with fields representing estimates)
    """)

    p.add_argument("infile", action="store", help="Input json")
    p.add_argument("--subset", action='store_true', help='Only run for representative subset of packages.')
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')
    p.add_argument("--method", default='bps', help='Estimation method to use')
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    p.add_argument("--ribologic_only", action='store_true', help='Run on RiboLogic FMN dataset only')
    p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_out.json')
    p.add_argument("--flanking", action='store_true', help='Run constructs with flanking sequences')
    args = p.parse_args()
    print(args)

    if args.parallel:
        pandarallel.initialize()

    df = pd.read_json(args.infile)

    print(df.keys())
    
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = args.infile.replace('.json','')

    if args.subset:
        if args.method=='Z':
            pkg_list = {k: packages_all_options[k] for k in ['vienna_2', 'rnastructure', 'contrafold_2', 'eternafold']}
        else:
            pkg_list = {k: packages_all_options[k] for k in subset_list}

    else:
        if args.method=='Z':
            pkg_list = packages_Z_subset
        else:
            pkg_list = packages_all_options


    if args.test:
        df = df.iloc[:1]

    if args.ribologic_only:
        df = df.loc[df['Round']=='Ribologic'][df['ligand']=='FMN']

    if args.method =='bps':

        control_values = {}
        for pkg_option, dct in pkg_list.items():
            control_values[pkg_option] = utils.compute_bp_kfold_ref_with_flanking(**dct)
            if args.verbose: print('wrote ref_kfold_est for %s' % pkg_option)

        for pkg_option, dct in pkg_list.items():
            if args.verbose: print('writing kfold_est_bp for %s' % pkg_option)
            #note: this is the inverse of kfold_est! Done to avoid handling the nan case if bp = 0 until later.

            if args.parallel and pkg_option != 'vienna_1':
                df['kfold_est_bp_%s' % pkg_option] = df.parallel_apply(lambda row: utils.compute_bp_kfold(row, flanking=args.flanking,**dct)/control_values[pkg_option], axis=1)
            else:
                df['kfold_est_bp_%s' % pkg_option] = df.apply(lambda row: utils.compute_bp_kfold(row, flanking=args.flanking, **dct)/control_values[pkg_option], axis=1)

            df.to_json("%s_est_bps_%s.json" % (outfile, todaysdate))

    elif args.method == 'Z':
        control_values = {}
        for pkg_option, dct in pkg_list.items():
            control_values[pkg_option] = utils.compute_Z_kfold_ref_with_flanking(**dct)
            if args.verbose: print('wrote ref_kfold_est for %s' % pkg_option)

        for pkg_option, dct in pkg_list.items():

            #note: produces 'kfold_est_Z' as a dict of two values: 'lig' and 'nolig'. Done to avoid computing Z's twice for sequences.
            if args.parallel and pkg_option != 'vienna_1':
                df[['log_kfold_unnorm_lig_Z_%s' % pkg_option,'log_kfold_unnorm_nolig_Z_%s' % pkg_option]] =\
                df.parallel_apply(lambda row: utils.compute_Z_kfold(row, flanking=args.flanking, **dct), axis=1,result_type="expand")
            else:
                df[['log_kfold_unnorm_lig_Z_%s' % pkg_option,'log_kfold_unnorm_nolig_Z_%s' % pkg_option]] =\
                df.apply(lambda row: utils.compute_Z_kfold(row, flanking=args.flanking, **dct), axis=1,result_type="expand")

            df['log_kfold_est_nolig_Z_%s' % pkg_option] =  df['log_kfold_unnorm_nolig_Z_%s' % pkg_option] - control_values[pkg_option]
            df['log_kfold_est_lig_Z_%s' % pkg_option] = df['log_kfold_unnorm_lig_Z_%s' % pkg_option] - control_values[pkg_option]

            df['log_AR_est_%s'% pkg_option] = np.where(df['switch']=='OFF',
             df['log_kfold_unnorm_lig_Z_%s' % pkg_option]-df['log_kfold_unnorm_nolig_Z_%s' % pkg_option],
              df['log_kfold_unnorm_nolig_Z_%s' % pkg_option]-df['log_kfold_unnorm_lig_Z_%s' % pkg_option])

            if args.verbose: print('wrote kfold Z estimations for %s' % pkg_option)
            df.to_json("%s_est_Z_%s.json" % (outfile, todaysdate))


