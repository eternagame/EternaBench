import numpy as np
import pandas as pd
import os, argparse
from arnie.bpps import bpps
from pandarallel import pandarallel
import RDAT_utils as utils
from datetime import date

todaysdate = date.today().strftime("%d%b%Y")

ARNIEDIR = '/home/users/hannahw1/arnie'

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

'vienna_langdon_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_langdon2018.par' % ARNIEDIR},
'vienna_rnasoft_pars': {'package': 'vienna_2', 'param_file':'%s/parameter_files/rna_andronescu2007.par' % ARNIEDIR},
'rnasoft_99' : {'package':'rnasoft_99'},
'rnasoft_07' : {'package':'rnasoft_07'},
'rnasoft_blstar' : {'package':'rnasoft'},
'rnasoft_99_nodangles' : {'package':'rnasoft_99-no-dangles'},
'rnasoft_bl_nodangles' : {'package':'rnasoft_bl-no-dangles'},
'rnasoft_lam-cg' : {'package':'rnasoft_lam-cg'},
'rnasoft_nom-cg' : {'package':'rnasoft_nom-cg'},

'eternafold_A': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_A' % ARNIEDIR},
'eternafold_B': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_B' % ARNIEDIR},
'eternafold_C': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_C' % ARNIEDIR},
'eternafold_D': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_D' % ARNIEDIR},
'eternafold_E': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_E' % ARNIEDIR},
'eternafold_F': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_F' % ARNIEDIR},
'eternafold_G': {'package': 'contrafold_2', 'param_file':'%s/parameter_files/eternafold_versions/params_G' % ARNIEDIR},

 }

most_recent_subset = ['vienna_2', 'nupack_99','rnastructure', 'contrafold_2', 'rnasoft_blstar','eternafold_B']

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Estimate reactivity as unpaired probability for sequences in input dataset.
    Input: .json file (written from rdats)
    Output: .json file (same as input, but with fields representing estimates)
    """)

    p.add_argument("infile", action="store", help="Input json")
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
    p.add_argument("--subset", action='store_true', help='Only run for one representative set of options for each package.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_out.json')

    args = p.parse_args()
    print(args)

    if args.parallel:
        pandarallel.initialize()

    df = pd.read_json(args.infile)

    if args.outfile:
        outfile = args.outfile
    else:
        outfile = os.path.basename(args.infile).replace('.json','')

    if args.test:
        df = df.iloc[:3]

    if args.subset:
        dct_to_run = {k:packages_all_options[k] for k in most_recent_subset}
    else:
        dct_to_run = packages_all_options

    for pkg_option, dct in dct_to_run.items():
        if args.verbose: print('writing p_unp for %s' % pkg_option)

        if args.parallel and pkg_option != 'vienna_1':
            df['full_p_%s' % pkg_option] = df.parallel_apply(lambda row: utils.write_unpaired_p(row, **dct), axis=1)
            df['p_%s' % pkg_option] = df.parallel_apply(lambda row: [row['full_p_%s' % pkg_option][x] for x in row['seqpos'] if x<len(row['sequence'])], axis=1)

        else:
            df['full_p_%s' % pkg_option] = df.apply(lambda row: utils.write_unpaired_p(row, **dct), axis=1)
            df['p_%s' % pkg_option] = df.apply(lambda row: [row['full_p_%s' % pkg_option][x] for x in row['seqpos'] if x<len(row['sequence'])], axis=1)

        df.to_json("%s_full_output_%s.json.zip" % (outfile, todaysdate))
