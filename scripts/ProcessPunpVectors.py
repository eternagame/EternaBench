import pandas as pd

from glob import glob
import argparse
import numpy as np
import sys
# from pandarallel import pandarallel
from tqdm import tqdm

tqdm.pandas()
# pandarallel.initialize(progress_bar=True)

def calculate_punp_variance(row, package_list):

    punp_arr = np.vstack([row[x] for x in package_list])
    # print(punp_arr.shape)
    # print(len(np.std(punp_arr, axis=0)))
    return np.median(np.std(punp_arr, axis=0))

    # ex_punp_vec = row[package_list[0]]
    # stdevs=[]
    # for i in range(len(ex_punp_vec)):
    #     tmp=[]
    #     for pkg in package_list:
    #         if len(row[pkg]) != len(ex_punp_vec):
    #             print('HERE', pkg, len(row[pkg]), len(ex_punp_vec))
    #         tmp.append(row[pkg][i])

    #     stdevs.append(np.std(tmp))

    # return np.mean(stdevs)


if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Process P(unp) vector files.
    Calculates variance over p(unp) values and re-splits by dataset.

    Input: keyword that all .json.zip's to compile should start with
    Output: .json file 
    """)

    p.add_argument("keyword", action="store", help="Keyword that all jsons to read in start with.")
    p.add_argument("--split_by", action="store", help="Field to re-split results on.")
    p.add_argument("-o", action='store', dest="outfile", help="Output name.")

    args = p.parse_args()

    files = sorted(glob("%s*.json.zip" % args.keyword))

    print('reading in file %s '% files[0])
    df = pd.read_json(files[0])

    for fil in files[1:]:
        print('reading in file %s '% fil)
        tmp = pd.read_json(fil)
        pkg_kws = [x for x in tmp.keys() if x.startswith('p_')]
        df = pd.concat([df, tmp[pkg_kws]],axis=1)

    package_list = [x for x in df.keys() if x.startswith('p_')]
    print("Packages found: ", package_list)
    df['punp_stdev'] = df.progress_apply(lambda row: calculate_punp_variance(row, package_list), axis=1)

    print(df.groupby('Dataset')['punp_stdev'].mean())
    print(df.groupby('project_name')['punp_stdev'].mean())

    if args.split_by:
        if args.split_by not in df.keys():
            raise RuntimeError("Data does not contain field %s" % args.split_by)
        print('Re-splitting results based on keyword %s' % args.split_by)
        uniq_split_keywords = list(df[args.split_by].unique())
        for kw in uniq_split_keywords:
            tmp = df.loc[df[args.split_by]==kw]
            tmp.to_json('%s_%s.json.zip' % (args.outfile, kw.replace(' ','_')))

    df = df.drop(columns=package_list)
    df.to_json('%s_punp_stdev.json.zip' % args.outfile)
            