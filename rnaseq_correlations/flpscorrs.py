from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import json

# custom modules
from bootstrap import *
from pairwise_corr import *

datadir = './input/'

# column names
names=['wbid','gname']
# load gene lists
allgenes = pd.read_csv(datadir + 'ce10genes.csv', header=None, names=names)
flps = pd.read_csv(datadir + 'flps.csv', header=None, names=names)
ins = pd.read_csv(datadir + 'ins.csv', header=None, names=names)
nlps = pd.read_csv(datadir + 'nlps.csv', header=None, names=names)


def bs_corr(gs, allgenes, exp_table, size=10000, save_to=False):
    """
    Boostrap test for correlation within a given gene set
    Generates (size) random replicates and computes the correlation within
    the set using exp table.

    gs: DataFrame, gene set of interest
    exp_table: str, path to DataFrame containing expression data
    size: int, number of bootstrap replicates
    """
    meancorr_r = np.empty(size)
    # keep only genes that appear in exp_table
    exp_table = pd.read_csv(exp_table)
    allgenes = allgenes[allgenes['wbid'].isin(exp_table['wbid'])]
    for i in range(size):
        print('bs iteration {}'.format(i))
        # generate random gene set
        rgs = allgenes.sample(len(gs)).reset_index(drop=1)
        # compute correlations within geneset
        corr_bs = pairwise_corr(rgs, rgs, exp_table)['spearmanr']
        # compute corrleations mean
        meancorr_r[i] = np.mean(corr_bs)

    # correlations of gene set
    corr_gs = pairwise_corr(gs, gs, exp_table)['spearmanr']
    # mean of gene set corrs
    meancorr_gs = np.mean(corr_gs)
    # p-value
    pval = np.sum(meancorr_r >= meancorr_gs) / len(meancorr_r)
    results = {'p-value':pval, 'mean_corr_gs':meancorr_gs,
                'mean_corr_random': list(meancorr_r)}
    if save_to:
        with open(save_to, 'w') as f:
            json.dump(results, f)
    return results

# using Dauer exp data
dauer = bs_corr(flps, allgenes,
            datadir + 'Dauerdata_meanTPM.csv',
            save_to='./output/dauer.txt')

# using Gerstein exp data
bs_gers, flps_gers, pval_gers = bs_corr(flps, allgenes,
                datadir + 'WBPaper00037953_expTPM_ce.csv', 
                save_to='./output/gerstein.txt')

flps_corrd = np.mean(pairwise_corr(flps, flps, datadir + 'Dauerdata_meanTPM.csv')['spearmanr'])
flps_corrg = np.mean(pairwise_corr(flps, flps, datadir + 'WBPaper00037953_expTPM_ce.csv')['spearmanr'])
