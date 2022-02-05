#!/usr/bin/env python

# tested with python v3.6.8
# guide scorer for cas12a capture

import math
import gzip
import numpy as np
import numpy.ma as ma
import scipy.stats
import pandas as pd
import argparse
from collections import Counter
from sklearn.linear_model import LinearRegression

# enumerate features for the oligo seq
def gen_feats(seq):
    
    base_pos_dict = {}
    for pos in list(range(0, 34)):
        base_pos_dict[str(pos)+'A'] = []; base_pos_dict[str(pos)+'C'] = []
        base_pos_dict[str(pos)+'G'] = []; base_pos_dict[str(pos)+'T'] = []
    for pos in range(0, len(seq)):
        base_pos_dict[str(pos)+seq[pos]].append(1)
        for alt_base in [i for i in ['A', 'G', 'T', 'C'] if i != seq[pos]]:
            base_pos_dict[str(pos)+alt_base].append(0)
    df = pd.DataFrame(base_pos_dict, index=[seq])
    GC_cont = (len([i for i in seq if i in ['G', 'C']]))/ len(seq)
    df['GC_cont'] = GC_cont
    
    GC_cont_only_guide = (len([i for i in seq[8:28] if i in ['G', 'C']])/ len(seq[8:28]))
    df['GC_cont_only_guide'] = GC_cont_only_guide
    
    num_GCs = (len([i for i in seq[8:28] if i in ['G', 'C']]))
    GC_imbalance = str(abs(np.subtract(0.5, (num_GCs/ 20))))
    df['GC_imbalance'] = GC_imbalance
    
    overhangGC = (len([i for i in seq[26:31] if i in ['G', 'C']]) / len(seq[18:23]))
    df['overhangGC'] = overhangGC
    
    dinuc_pos_dict = {}
    dinuc_list = ['AA', 'AC', 'AG', 'AT', 'CC', 'CA', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    
    for pos in list(range(0, 34-1)):
        for dinuc in dinuc_list:
            dinuc_pos_dict[str(pos)+dinuc] = []
    
    for pos in range(0, 34-1):
        dinuc_pos_dict[str(pos)+seq[pos:pos+2]].append(1)
        for alt_dinuc in [i for i in dinuc_list if i != seq[pos:pos+2]]:
            dinuc_pos_dict[str(pos)+alt_dinuc].append(0)
    df = pd.concat([df, pd.DataFrame(dinuc_pos_dict, index = df.index)], axis=1)
    return df

def main():

    p = argparse.ArgumentParser(description='Score guides for cas12a capture. Returns score to standard output')

    p.add_argument('-s', '--seq', action='store', required=True, 
        type=str, metavar='ACGT',
        help='list of input sequences, one per line')

    p.add_argument('-f', '--features', action='store', required=True, 
        type=str, metavar='feature_file',
        help='feature file')

    p.add_argument('-t', '--training', action='store', required=True, 
        type=str, metavar='training_data',
        help='training data file')

    args = p.parse_args()

    # read sequences
    seqs = []
    with open (args.seq) as sequences:
       for line in sequences:
           ls = line.strip()
           seqs.append(ls)

    # for each sequence
    for seq in seqs:

        # select feature
        query_df = gen_feats(seq)

        selected_feats = []
        with open (args.features) as data:
            for line in data:
                ls = line.strip()
                selected_feats.append(ls)

        query_df_selected_feats = query_df[selected_feats]

        # get training data
        df = pd.read_csv(args.training, index_col = 0)

        # add a pseudocount to all guides so that we can still take log:
        df['sum'] = df['sum']+1

        # calculate log sum:
        df['log_sum'] = np.log10(df['sum'])

        # train the model
        X_train = df[selected_feats].copy()
        y_train = df['log_sum']
        regr = LinearRegression().fit(X_train, y_train)

        # make predictions
        preds = regr.predict(query_df_selected_feats)
        print(seq,"\t",preds)

if __name__=='__main__':
  main()
