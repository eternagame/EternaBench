import numpy as np
import pandas as pd
import os, argparse
import arnie
from arnie.utils import load_package_locations
from arnie.bpps import bpps
from pandarallel import pandarallel
import eternabench.chemmapping_utils as utils
from eternabench.package_utils import package_options, example_subset

from datetime import date

todaysdate = date.today().strftime("%d%b%Y")

def write_unpaired_p(row, package='vienna_2', DEBUG=False, **kwargs):
    try:
        vec = 1 - np.sum(bpps(row['sequence'], package=package, DEBUG=DEBUG, **kwargs), axis=0)
    except:
        vec = np.NaN*np.zeros(len(row['sequence']))
    return vec

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Estimate reactivity as unpaired probability for sequences in input dataset.
    Input: .json file (written from rdats)
    Output: .json file (same as input, but with fields representing estimates)
    """)

    p.add_argument("infile", action="store", help="Input json")
    p.add_argument("--package", action='store', help='package keyword in Arnie (see doc for full list of options).')
    p.add_argument("--test", action='store_true', help='Tests first 3 constructs in each package.')
    p.add_argument("--subset", action='store_true', help='(DEP) Run for one representative set of options for each package.')
    p.add_argument("--parallel", action='store_true', help='Runs in parallel using Pandarallel.')
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")
    p.add_argument('--debug', action="store_true", help='Print full output from arnie.')
    p.add_argument("-o", action="store", dest='outfile', help='name of output json file (default is <input_name>_out.json')

    args = p.parse_args()

    print("#### ETERNABENCH: Calculate P(Unpaired) Vectors #####")
    print("Last Update July 2021, H Wayment-Steele")
    print("Running calculations using Arnie located at %s" % arnie.__path__)
    print("Arnie paths:")
    arnie_dct = load_package_locations()
    for k, v in arnie_dct.items():
        print("%s: %s" % (k,v))

    if args.parallel:
        pandarallel.initialize()

    df = pd.read_json(args.infile)

    df = df.dropna(subset=['sequence'])
    df['sequence'] = [x.replace('.','').replace('\n','') for x in df['sequence']]

    if args.outfile:
        outfile = args.outfile
    else:
        outfile = os.path.basename(args.infile).replace('.json.zip','')+'_'

    if args.test:
        df = df.iloc[:3]

    if args.test:
        packages_to_run = package_options
    else:
        packages_to_run = [args.package]

    for pkg_option in packages_to_run:

        if args.verbose: print('writing p_unp for %s' % pkg_option)

        if args.parallel and pkg_option != 'vienna_1':
            if 'full_p_%s' % pkg_option not in df.keys():
                df['full_p_%s' % pkg_option] = df.parallel_apply(lambda row: write_unpaired_p(row, DEBUG=args.debug, **package_options[pkg_option]), axis=1)
            else:
                print('Found %s' % pkg_option)
            df['p_%s' % pkg_option] = df.parallel_apply(lambda row: [row['full_p_%s' % pkg_option][x] for x in row['seqpos'] if x<len(row['sequence'])], axis=1)
            df = df.drop(columns=['full_p_%s' % pkg_option])

        else:
            if 'full_p_%s' % pkg_option not in df.keys():
                df['full_p_%s' % pkg_option] = df.apply(lambda row: write_unpaired_p(row, DEBUG=args.debug, **package_options[pkg_option]), axis=1)
            else:
                print('Found %s' % pkg_option)
            df['p_%s' % pkg_option] = df.apply(lambda row: [row['full_p_%s' % pkg_option][x] for x in row['seqpos'] if x<len(row['sequence'])], axis=1)
            df = df.drop(columns=['full_p_%s' % pkg_option])

        if not args.test: 
            df.to_json("%s_%s.json.zip" % (outfile, args.package))
        else:
            df.to_json("%s.json.zip" % (outfile))


    # outdir='%s_Dataset_chunks_%s' % (outfile, todaysdate)

    # print('Writing output, chunked by Dataset, to %s' % outdir)

    # if not os.path.exists(outdir):
    #     os.mkdir(outdir)

    # for filename in df.Dataset.unique():
    #     tmp = df.loc[df.Dataset==filename]
    #     tmp.to_json('%s/%s.json.zip' % (outdir, filename.replace(' ','_')))


