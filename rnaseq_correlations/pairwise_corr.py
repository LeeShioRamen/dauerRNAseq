"""
Compute gene pairwise correlation in RNAseq expression.
Takes two file names as arguments from stdin, which must be comma separated
value files without header, with rows of: WBID, name.

Author: Porfirio Quintero-Cadena
"""

import sys
import os
import pandas as pd

def add_exp_pairs(df, exp_table, ID, FPKM_to_TPM): 
    """
    get multiple expression indices for left and right genes of pairs

    Arguments
    -------

    df: dataframe with paired genes
    exp_table: path to csv file with expression data in FPKM or TPM
    ID: name of gene_ID column (according to organism, e.g. wbid, FBgn, ENSG)
    FPKM_to_TPM: normalize each RNAseq exp to TPM as: FPKM_gene/FPKM_sum*1e6

    Returns
    -------
    gene_1, gene_2 : Dataframes with gene_ID, exp data for l and r gene in pair


    """
    # Load table of gene expression and WBID
    if type(exp_table) is str:
        exp_table = pd.read_csv(exp_table)
    exp_table = exp_table.apply(pd.to_numeric, errors="ignore")
    exp_numeric = exp_table._get_numeric_data()
    numeric_cols = exp_numeric.columns

    # convert to TPM from FPKM if requested
    if FPKM_to_TPM:
        exp_table[numeric_cols] = exp_numeric/exp_numeric.sum(axis=0)*1e6

    # delete genes that are not detected in any experiment
    # or detected in less than 80% of experiments
    never_expressed = exp_table[(exp_numeric.sum(axis=1)==0)|
        (exp_numeric.count(axis=1)<exp_numeric.shape[1]*0.8)][ID]
    exp_table = exp_table[~(exp_table[ID].isin(never_expressed))]

    # transformed expression values to ranked values for spearman correlation
    exp_table[numeric_cols] = exp_table._get_numeric_data().rank(axis=0)
    gene_1 = pd.DataFrame()
    gene_2 = pd.DataFrame()
    gene_1[ID] = df[ID+'1']
    gene_2[ID] = df[ID+'2']
    gene_1 = pd.merge(gene_1, exp_table, how='left', on=ID)
    gene_2 = pd.merge(gene_2, exp_table, how='left', on=ID)

    return gene_1, gene_2

def compute_corr(self, other):
    """
    Compute pairwise spearman correlation between rows of
    two DataFrame objects.

    Arguments
    ----------
    other : DataFrame

    Returns
    -------
    correls : Series
    """
    # compute pearson correlation 
    # This would be spearmanr if values were ranked as using add_exp_pairs
    this = self._get_numeric_data()
    other = other._get_numeric_data()

    return this.corrwith(other, axis=1) 

def pairwise_corr(gene_list1, gene_list2, exp_dir, remove_self_pairs=True):
    """
    Computes pairwise correlation between two lists of genes

    Arguments
    ----------
    gene_list1, gene_list2: lists of genes, each gene with WBID and name
    remove_self_pairs: wheter to remove pairs of the same gene

    Returns
    -------
    pairs : DataFrame with gene pairs and spearman correlation

    """
    genes = gene_list1.values.tolist()
    pairs = pd.DataFrame(columns=['wbid1', 'gname1','wbid2', 'gname2'])
    # generate all pairs
    for gene in genes:
        curr_pair = pd.DataFrame(columns=['wbid1', 'gname1','wbid2', 'gname2'])
        curr_pair[['wbid2','gname2']] = gene_list2
        curr_pair[['wbid1','gname1']] = gene
        pairs = pairs.append(curr_pair)
    pairs.index = range(len(pairs))

    print('adding expression data...')
    g1_exp, g2_exp = add_exp_pairs(pairs, exp_table=exp_dir, ID='wbid', 
            FPKM_to_TPM=False)
    print('computing correlation...')
    spearmanr = compute_corr(g1_exp, g2_exp)
    pairs['spearmanr'] = spearmanr
    print('done')
    if remove_self_pairs: pairs = pairs[(pairs.wbid2 != pairs.wbid1)]
    return pairs
