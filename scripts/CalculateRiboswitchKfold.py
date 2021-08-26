import numpy as np
import pandas as pd
import os, argparse
import arnie
from arnie.utils import load_package_locations

from arnie.bpps import bpps
from pandarallel import pandarallel
import eternabench.riboswitch_utils as utils
from datetime import date
from eternabench.package_utils import package_options, example_subset
from tqdm import tqdm

tqdm.pandas()

todaysdate = date.today().strftime("%d%b%Y")

packages_Z_subset = {k: package_options[k] for k in package_options.keys()
     if k.startswith('vienna') or k.startswith('contrafold') or k.startswith('eternafold') or k.startswith('rnastructure')}

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Estimate protein affinity for sequences in input dataset.
    Input: .json file (eterna or ribologic dataset)
    Output: .json file (same as input, but with fields representing estimates)
    """)

    p.add_argument("infile", action="store", help="Input json")
    p.add_argument("--package", action='store', help='package keyword in Arnie (see doc for full list of options).')
    #p.add_argument("--subset", action='store_true', help='Only run for representative subset of packages.')
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')
    p.add_argument("--method", default='bps', help='Estimation method to use')
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    p.add_argument("--ribologic_only", action='store_true', help='Run on RiboLogic FMN dataset only')
    p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_out.json')
    p.add_argument("--flanking", action='store_true', help='Run constructs with flanking sequences')
    args = p.parse_args()
    print(args)

    print("#### ETERNABENCH: Calculate Riboswitch Kfold #####")
    print("Last Update July 2021, H Wayment-Steele")
    print("Running calculations using Arnie located at %s" % arnie.__path__)
    print("Arnie paths:")

    arnie_dct = load_package_locations()
    for k, v in arnie_dct.items():
        print("%s: %s" % (k,v))

    if args.parallel:
        pandarallel.initialize(progress_bar=True)

    df = pd.read_json(args.infile)

    print(df.keys())
    
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = args.infile.replace('.json','')

    # if args.subset:
    #     if args.method=='Z':
    #         pkg_list = {k: package_options[k] for k in ['vienna_2', 'rnastructure', 'contrafold_2', 'eternafold']}
    #     else:
    #         pkg_list = {k: package_options[k] for k in subset_list}

    # else:
    #     if args.method=='Z':
    #         pkg_list = packages_Z_subset
    #     else:
    #         pkg_list = package_options


    if args.test:
        df = df.iloc[:3]
        packages_to_run = package_options
    else:
        packages_to_run = [args.package]

    if args.ribologic_only:
        df = df.loc[df['Round']=='Ribologic'][df['ligand']=='FMN']

    if args.method =='bps':

        control_values = {}
        for pkg_option in packages_to_run:
            control_values[pkg_option] = utils.compute_bp_kfold_ref_with_flanking(**package_options[pkg_option])
            if args.verbose: print('wrote ref_kfold_est for %s' % pkg_option)

        for pkg_option in packages_to_run:
            if 'kfold_est_bp_%s' % pkg_option not in df.keys():
                if args.verbose: print('writing kfold_est_bp for %s' % pkg_option)
                if args.parallel and pkg_option != 'vienna_1':
                    # log (p_ref / p_ij)
                    df['log_kfold_est_bp_%s' % pkg_option] = df.parallel_apply(lambda row: np.log(control_values[pkg_option]/utils.compute_bp_kfold(row, flanking=args.flanking,**package_options[pkg_option])), axis=1)
                else:
                    df['log_kfold_est_bp_%s' % pkg_option] = df.progress_apply(lambda row: np.log(control_values[pkg_option]/utils.compute_bp_kfold(row, flanking=args.flanking, **package_options[pkg_option])), axis=1)
            else:
                print('found %s' % pkg_option)

        if not args.test: 
            df.to_json("%s_%s_bps.json.zip" % (outfile, args.package))
        else:
            df.to_json("%s_bps.json.zip" % (outfile))

    elif args.method == 'Z':
        control_values = {}
        for pkg_option in packages_to_run:
            control_values[pkg_option] = utils.compute_Z_kfold_ref_with_flanking(**package_options[pkg_option])
            if args.verbose: print('wrote ref_kfold_est for %s' % pkg_option)

        for pkg_option in packages_to_run:
            if 'log_kfold_unnorm_lig_Z_%s' % pkg_option not in df.keys():
                #note: produces 'kfold_est_Z' as a dict of two values: 'lig' and 'nolig'. Done to avoid computing Z's twice for sequences.
                if args.parallel and pkg_option != 'vienna_1':
                    df[['log_kfold_unnorm_lig_Z_%s' % pkg_option,'log_kfold_unnorm_nolig_Z_%s' % pkg_option]] =\
                    df.parallel_apply(lambda row: utils.compute_Z_kfold(row, flanking=args.flanking, **package_options[pkg_option]), axis=1,result_type="expand")
                else:
                    df[['log_kfold_unnorm_lig_Z_%s' % pkg_option,'log_kfold_unnorm_nolig_Z_%s' % pkg_option]] =\
                    df.progress_apply(lambda row: utils.compute_Z_kfold(row, flanking=args.flanking, **package_options[pkg_option]), axis=1,result_type="expand")

                df['log_kfold_est_nolig_Z_%s' % pkg_option] =  df['log_kfold_unnorm_nolig_Z_%s' % pkg_option] - control_values[pkg_option]
                df['log_kfold_est_lig_Z_%s' % pkg_option] = df['log_kfold_unnorm_lig_Z_%s' % pkg_option] - control_values[pkg_option]

                df['log_AR_est_%s'% pkg_option] = np.where(df['switch']=='OFF',
                df['log_kfold_unnorm_lig_Z_%s' % pkg_option]-df['log_kfold_unnorm_nolig_Z_%s' % pkg_option],
                df['log_kfold_unnorm_nolig_Z_%s' % pkg_option]-df['log_kfold_unnorm_lig_Z_%s' % pkg_option])

                if args.verbose: print('wrote kfold Z estimations for %s' % pkg_option)
            else:
                print('found %s' % pkg_option)
        if not args.test: 
            df.to_json("%s_%s_Z.json.zip" % (outfile, args.package))
        else:
            df.to_json("%s_Z.json.zip" % (outfile))

