import pandas as pd
import argparse
from glob import glob
from eternabench.stats import calculate_Z_scores


import sys


if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    '''Compile results from EternaBench bootstrapping.
    Input: keyword that all .json's to compile should start with
    Output: .json file''')

    p.add_argument("keyword", action="store", help="Keyword that all jsons to read in start with.")
    p.add_argument("-o", action='store', dest="outfile", help="Output name.")
    p.add_argument('--metric', action='store', default='pearson', help='name of metric field, default=pearson.')
    p.add_argument('--dataset_field', action='store', default='Dataset', help="name of dataset field, default='Dataset'")
    p.add_argument('--ranking_category', action='store', default=None, help="If provided, splits ranking into these classes.")
    p.add_argument('--calculate_Z_scores', action='store', help="Calculate Z scores based on package list provided.")
    p.add_argument('--dataset_list', action='store', help='Only calculate over specified datasets.')
    args = p.parse_args()

    files = sorted(glob("%s*_BOOTSTRAPS.json.zip" % args.keyword))

    print("Found %d files" % len(files))

    df = pd.DataFrame()

    for fil in files:
        #print('reading in file %s '% fil)
        tmp = pd.read_json(fil)
        df = df.append(tmp, ignore_index=True)

    #df.to_json('%s.json.zip' % args.outfile)

    if args.dataset_list is not None:
        dataset_list = []
        with open(args.dataset_list,'r') as f:
            for lin in f.readlines():
                dataset_list.append(lin.strip())

        print("Calculating for datasets: ", dataset_list)

        df = df.loc[df[args.dataset_field].isin(dataset_list)]

    if args.calculate_Z_scores is not None:
        package_list=[]
        with open(args.calculate_Z_scores,'r') as f:
            for lin in f.readlines():
                package_list.append(lin.strip())


        print("Read in packages from package list: ", package_list)

        sorted_zscore_stats, ranking = calculate_Z_scores(df, package_list, metric=args.metric, dataset_field=args.dataset_field, ranking_category=args.ranking_category, include_efold=True)

        sorted_zscore_stats.to_csv('%s_%s_zscores_by_%s.csv' % (args.outfile, args.metric, args.dataset_field), index=False)
        ranking.to_csv('%s_%s_ranking.csv' % (args.outfile, args.metric), index=False)
